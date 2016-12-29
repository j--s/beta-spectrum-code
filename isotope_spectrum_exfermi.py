from scipy.integrate import quadrature
from scipy.integrate import romb
from scipy.integrate import quad
from scipy.special import gamma
from operator import add
import numpy as np
import time
import sys

start_time = time.time()
procfile = open('proc_err.dat', 'w+')
def read_prod_data(file):

    """Reads the ENDF product yield data file"""

    prod_file = open(file,'r')      
    
    pfile = prod_file.readlines()
    Z = []
    A = []
    isomer = []
    element = []
    prod_yield = []
    for lines in pfile:
        word = lines.split() 
        Z.append(float(word[1]))
        element.append(word[2])
        A.append(float(word[3]))
        isomer.append(int(float(word[5])))
        prod_yield.append(float(word[6]))

    return A, Z, element, isomer, prod_yield 

def read_decay_data(decay_file, parent, isomer_flag):

    """Reads the extracted decay data from the betadata.dat file."""

    print('------------------------------------------', file=procfile)
    print('Parent: ', parent, file=procfile)
    isomer = 0
    parent_e = []
    decay_data = []
    level_decay_data = []

    f = open(decay_file, 'r')
    lines = f.readlines()
    f.close()

    i = 0
    ii = 1
    for line in lines:
        objects = line.split()

        if  objects[0] == parent:
            if objects[1] == isomer_flag: 
                print('Isomer Flag: ', isomer_flag, file=procfile)
                i = 1
        elif objects[0] == 'EOF' and (i == 1 or i == 2):
            decay_data = level_decay_data
            return parent_e, decay_data
        elif i == 1:
            print('Parent E: ', objects[0], file=procfile)
            parent_e = [float(objects[0]), float(objects[1])]
            i = 2
        elif i == 2:
            level_decay_data.append(objects)
        elif ii > len(lines) - 1:
            decay_data = level_decay_data
            return parent_e, decay_data

        ii += 1

    return parent_e, decay_data

def fermi(a, z, e, e0, forbid):

    """Calculates the Fermi function for a beta- decay"""

    alpha = 0.007297                    #Fine structure constant from NIST website  
    rad =  0.0031 * (a**(1/3))             #Nuclear charge radius
    p = np.sqrt(e*e + 2 * 511. * e)    #Electron relativistic energy
    y = (alpha * z * (e + 511.)) / p 
    gam = np.sqrt(1. - (alpha * z)**2)
    num = gamma(gam + y * 1j)           #Calculating the gamma functions
    denom = gamma(2. * gam + 1.)
    frac = np.square(np.absolute(num) / np.absolute(denom))

    c = correction(e, e0, forbid)
    
    #Combine everything into the Fermi function     
    L_0 = 0.5 * (1 + gam)

    F_0 = 4. *  ((2. * p * rad)**(-2. * (1. - gam))) * frac * c * np.exp(np.pi * y)                                   
    f = F_0 * L_0                                                                                                     

    return f 

def correction(e, e0, forbid):

    """Shape factors for allowed (A), first-forbidden (FF), and 
    first-forbidden unique (FFU) decays"""

    #Defining the electron and neutrino momentum
    p = np.sqrt(e*e + 2 * 511. * e) 
    q = (e0 - e)     
    beta = p / e    
    #The shape factors, determined by the spin-parity of decay
    #Gamow-Teller decays
    if forbid == 'A':
        c = 1
    elif forbid == 'FF': 
        c = p*p + q*q + 2 * beta*beta * q * e
    elif forbid == 'FFU':
        c = p*p + q*q 
    #Find an expression for this forbideness
    elif forbid == 'SFU':
        c = 1
    #    #c = p*p - 2 * e * q
    #Fermi decays
    #elif forbid == '0+':  
    #    c = 1
    #elif forbid == '1-':
    #    c = p*p + q*q + (2/3) * beta*beta * q * e

    return c

def probability(a, z, e, e0, forbid):

    """Function to calculate the unormalized probability spectrum 
    of beta- decay branch. Includes only the finite size correction... so far."""

    alpha = 0.007297
    p = np.sqrt(e*e + 2 * 511. * e)
    q = (e0 - e)**2
    rad = 0.0031 * (a**(1/3))
    #Finite size correction
    delta = ((- 1.6) * z * alpha * rad * e) * (1 + ((0.3214 * np.square(511.)) / (e*e)))
    fermi_func = fermi(a, z, e, e0, forbid)
    #Calculate the beta spectrum.
    prob = (p * q * (e + 511.) * fermi_func * (1 + delta))

    return prob

def normalization(a, z, e0, forbid):
    
    """Integrates the analytical form of the beta spectrum and returns the 
    normalization factor."""

    #Check out Scipy Ctype integration with C via Python, >2x speedup.
    y = lambda e: probability(a, z, e, e0, forbid)
    l_lim = 1.0
    norm = quad(y, l_lim, e0)[0]
    print('Branch Integration: ', norm, file=procfile)

    return norm

def total_normalization(spectrum):
    
    """Integrates the coordinate pairs of the total beta spectrum and 
    returns the normalization factor."""

    y1 = []
    for coord in spectrum:
        y1.append(coord[1])
    y = np.array(y1)
    #Check out Scipy Ctype integration with C via Python, >2x speedup.
    #Important!! The romb (romberg) integration must have 2^k + 1 points.
    l_lim = 1.0
    dx = 1
    norm = romb(y, dx=dx)

    return norm

def branch_spectrum(a, z, e0, forbidden, intensity):
    
    """Calculate the spectrum for a single decay branch. It includes the intensity factor."""

    branch_spectrum = []
    print('Forbiddeness: ', forbidden, file=procfile)
    I = intensity / 100
    print('Intensity: ', I, file=procfile)
    norm = normalization(a, z, e0, forbidden)
    #Create spectrum, note for the romb integration must have 2^k + 1 points.
    for i in range(0, 4097): 
        e = (i * 3.5) + 1.0#0.5
        if e <= e0:
            prob = probability(a, z, e, e0, forbidden)
            branch_spectrum.append([e, I * (prob / norm)])
        else: 
            branch_spectrum.append([e, 0])

    return branch_spectrum

def decay_spectrum(f, spectra, intensity):

    """Somewhat generic function to sum spectra, write to file,
    and return the spectrum sum.

    f = 0 only sums the spectra
    f = file sums and prints the complete spectrum to file."""

    #Sum spectra
    i = 1
    decay_spectrum = []
    for spectrum in spectra:
        print(len(spectrum), i, file=procfile)
        spec = []
        if i < 2:
            for prob in spectrum:
                decay_spectrum.append(prob[1])
        else:
            for total, prob in zip(decay_spectrum, spectrum):
                spec.append(total + prob[1])
            decay_spectrum = spec
        i += 1

    spectrum = spectra[0]
    isotope_spectrum = []

    #Multiply by intensity only
    if not f:
        for e, y in zip(spectra[0], decay_spectrum):     
            isotope_spectrum.append([e[0], y*intensity])
    #Multiply by ntensity and write to file
    else:
        for e, y in zip(spectra[0], decay_spectrum):     
            print('%2.4f   %1.20f ' % (e[0], y*intensity), file=f)
            isotope_spectrum.append([e[0], y*intensity])
      
    return isotope_spectrum


#Assign and open files to be read and written.
fission_file = 'cum_py.txt'
decay_file = 'betadata.dat'
h = open('final_spectrum.dat', 'w+')

#Gather the fission yield data
A, Z, element, isomer_flags, prod_yield = read_prod_data(fission_file)

ii = 0
isotope_spectrum = []
for a, z, elem, isomer_flag, pyield in zip(A, Z, element, isomer_flags, prod_yield):

    print('Fission Intensity: ', pyield, file=procfile) 
    parent = str(int(a)) + elem.upper()
    print(ii, ' ', parent)
    #Gather decay data
    parent_energy, level_decay_data = read_decay_data(decay_file, parent, str(isomer_flag))        
    level_spectrum = []                                                                            
    #Each ground state or isomer                                                                         
    print('Number of Branches: ', len(level_decay_data), file=procfile)                                           
    if len(level_decay_data) < 1:                                                                  
        print('No decay Branches!', len(level_decay_data), file=procfile)                                   
        ii += 1                                                                                    
        continue                                                                                   
    else:                                                                                          
        print(ii, ' ','Spectrum processed', parent)                                                    

        #file = 'beta_spec_test' + parent + str(ii) + '.dat' 
        #f = open(file, 'w+')                                                                       
        #Make list containing each branch spectrum.                                                
        spectra = []                                                                               
        #Calculate the spectrum for each branch                                                    
        for branch in level_decay_data:                                                            
            e_level = float(branch[0])                                                             
            forbidden = branch[1]                                                                  
            intensity = float(branch[2])                                                           
            #Disregard branches with intensity < 0.25%                                                
            if intensity < 0.25:                                                                   
                continue                                                                           
            else:                                                                                  
                e0 = parent_energy[1] + parent_energy[0]                                           
                e_max = (e0 - e_level)# / 1000  #Max energy of the electron in keV                 
                if e_max <= 0:                                                                     
                    continue                                                                       
                else:                                                                              
                    print('Daughter Energy: ', e_max, file=procfile)                                                     
                    spectra.append(branch_spectrum(a, z, e_max, forbidden, intensity))             
                                                                                                   
        print('Len of spectra', len(spectra), file=procfile)                                        
        if len(spectra) > 0:                                                                       
            #Sum and normalize the branch spectra into the level spectrum                          
            spectrum = decay_spectrum(0, spectra, 1.0)                                             
            level_norm = total_normalization(spectrum)                                             
            print('Nucleus or Isomer Integration: ', level_norm, file=procfile)                               
            level_inten = 1 / level_norm                                                           
            level_spectrum.append(decay_spectrum(0, spectra, 1.0))                                 
                                                                                                   
        #f.close()                                                                                  
                                                                                                   
    #Multiply the level spectrum by the fission yield factor (pyield).
    level_spectrum2 = []                                                                           
    if len(level_spectrum) > 0:                                                                    
        print('Multiplied by Fission Intensity', file=procfile)
        level_spectrum2 = (decay_spectrum(0, level_spectrum, pyield))                              
        isotope_spectrum.append(level_spectrum2)                                                   
                                                                                                
    ii += 1                                                                                        

#Sum and normalize all isotope spectra into final spectrum and write to final_spectrum_myfermi.dat file
fission_spectrum = []
fission_spectrum.append(decay_spectrum(0, isotope_spectrum, 1.0))
fission_norm = total_normalization(fission_spectrum)
print('Final Spectrum Integration: ', fission_norm, file=procfile)
intensity = 1 / fission_norm
decay_spectrum(h, isotope_spectrum, intensity) 

print('Number of files: ', ii)
h.close()
procfile.close()
print('Total run time:  ', time.time() - start_time)
