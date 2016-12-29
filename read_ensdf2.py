import numpy as np
import sys

test_file = open('read_test.dat','w+')

def read_prod_data(file):

    """Reads the ENDF product yield data file"""
    prod_file = open(file,'r')
   
    pfile = prod_file.readlines()
    prod_file.close()

    A = []
    element = []
    for lines in pfile:
        word = lines.split()
        element.append(word[2])
        A.append(word[3])

    return A, element

def read_yield_data(fission_file):

    """Reads the decay data extracted from the ENSDF file"""
    f = open(fission_file, 'r')
    lines = f.readlines()
    f.close

    intensity= []
    element = []
    a = []
    z = []
    for line in lines:
        elements.append(lines.split())
        a.append(elements[3])
        element.append(elements[2])
        z.append(elements[1])
        intensity.append(elements[6])

    return a, z, element, intensity

def data(file, errfile):
    
    """Function to extract Q_0, energy levels, total spin, and parity
    for each beta decay branch. It is designed to read from ENDSF
    formatted files."""

 
    print('------------------------------------------', file=errfile)
    print('------------------------------------------', file=errfile)

    #Read data from the ENSDF file 
    try:
        f = open(file, 'r')
    except:
        print('No file! ', file, file=errfile)
        return
    length = len(f.readlines())
    print('Length of file: ', length)
    f = open(file, 'r')
    e = []
    e2 = []
    inten = []
    forbidden = []
    norm = 1.
    i = 1
    ii = 10
    jj = 0
    kk = 0
    isomer = 0
    parent = ''
    energy_flag = 0

    while ii > 0:

        line = f.readline()
        ii = len(line.split())
        word = line.split()                      
        word2 = str(line)

        if not word: #End of isotope beta data
            total_inten = 0
            inten_norm = 1
            missed_intensity = 0
            missed_flag = 0
            print('Len of inten', len(inten))
            if len(inten) > 1:
                print('-----------', file=errfile)
                print('Branching Ratios', file=errfile)
                for intensities in inten:
                    print(intensities, file=errfile)
                    total_inten += intensities
                    if intensities == 0.:
                        missed_flag = 1
                        missed_intensity += 1
                print('Total Branching Ratio: ', '%3.2f' % total_inten, file=errfile)
                if total_inten > 101.:
                    print('Error!', '%3.2f' % total_inten, file=errfile)
                    print('Total Branching Ratio Greater than 100%', file=errfile)
                print('-----------', file=errfile)
                inten_residue = 100 - total_inten
                if inten_residue > 0 and missed_flag:
                    inten_norm = inten_residue / missed_intensity

            for e1, forbid1, inten1 in zip(e, forbidden, inten):
                if inten1 == 0. and missed_intensity:
                    print('%5.4f' % e1, '  ', forbid1, '    %2.4f  ' % (inten_norm), file=file2)
                else:
                    print('%5.4f' % e1, '  ', forbid1, '    %2.4f  ' % (inten1), file=file2)
                
            if i < length:
                isomer += 1
                e = []                             
                inten = []
                forbidden = []
                norm = 1.
                print('EOF', file=file2)
                ii = 10
        elif i < 2:
            daughter = word[0]
        else:
            if len(word) > 1:
                identifier = word[1]         
            else:
                identifier = ''
            #Extract normalization multiplier
            
            if identifier == 'N':           
                if word2[41:49] == '        ':
                    norm = 1.0
                else:
                    norm = float(word2[41:49])
                print('File: ', word[0], file=errfile)   
                                                                                    
            #Extract energy of parent nuclei
            elif identifier == 'P':
                
                if (i - jj) <= 1:
                    i += 1
                    continue
                elif  word2[64:74].split()[0] == '':
                    e0 = 0.
                    parent = word[0]
                    print('Parent: ', word[0], file=read_test)
                    print(parent, '  ', isomer, file=file2)
                else:                         
                    e0 = float(word2[64:74])
                    parent = word[0]
                    print('Parent: ', word[0])
                    print(parent, '  ', isomer, file=file2)
                try:                                       
                    parent_e = float(word2[9:19])
                except:                                    
                    parent_e = 0.
                print('%5.4f' % parent_e, '%5.4f ' % e0, file=file2)
                jj = i

            elif word[0] != daughter:                                                                                     
                i += 1
                continue

            #Extract energy level
            elif identifier == 'L':                          
                energy_flag = 1
                if (i - kk) <= 1:                
                    i += 1
                    continue
                elif len(e) > (len(inten) or len(forbidden)):       
                    e.pop()                                    
                    try:                                       
                        #Some energy levels have the form E+X or X+E, this mess rids the energy of '+X'
                        if str(word2[9:19].split()[0])[-1] == 'X':                                                                     
                            e.append(float(str(word2[9:19].split()[0])[0:(len(str(word2[9:19].split()[0]))-2)]))
                        elif str(word2[9:19].split()[0])[0] == 'X':
                            e.append(float(str(word2[9:19].split()[0])[2:(len(str(word2[9:19].split()[0])))]))
                        else:
                            e.append(float(word2[9:19]))
                    except:                                    
                        print('No level energy on line ', i, file=errfile) 
                else:                                          
                    try:                                       
                        if str(word2[9:19].split()[0])[-1] == 'X':
                            e.append(float(str(word2[9:19].split()[0])[0:(len(str(word2[9:19].split()[0]))-2)]))
                        elif str(word2[9:19].split()[0])[0] == 'X':
                            e.append(float(str(word2[9:19].split()[0])[2:(len(str(word2[9:19].split()[0])))]))
                        else:
                            e.append(float(word2[9:19]))
                    except:  
                        print('No level energy on line ', i, file=errfile)
                kk = i

            #Extract energy level, intensity, and decay variety for daughter nuclei
            elif identifier == 'B' and energy_flag:
                energy_flag = 0
                if len(e) <= (len(inten) or len(forbidden)):       
                    i += 1
                    continue
                else:
                    try:                                                                                                 
                        e2.append(float(word2[9:19]))
                    except:                                    
                        e2.append(0.)
                    try:                                                           
                        if word2[21:28].split()[0][0] == '(':
                            inten.append(float(word[21:28].split()[0][1:(len(word[21:28].split()[0])-1)])*norm)
                        else:
                            inten.append(float(word2[21:28])*norm)                     
                    except:                                                        
                        print('No Branching Ratio on line: ', i, file=errfile)
                        inten.append(0.)
                                                                                                                
                    if word[-1] == '?':                                               
                        print('Line', (i), 'Uncertain Beta decay',file=errfile)
                        forbidden.append('A')
                    elif word[-1] == 'S':                                            
                        print('Line', (i), 'Predicted Beta decay',file=errfile)
                        forbidden.append('A')
                    elif word[-1] == '1U':                                           
                        forbidden.append('FFU')                                    
                    elif word[-1] == '2U':
                        print(i, '  ', word[-1])
                        forbidden.append('SFU')                                    
                    elif word[-1] == '1':
                        forbidden.append('FF')                                     
                    elif word2[77:79] == '  ':
                        forbidden.append('A')                                      
                    
                    
        i += 1

    i -= 1


    if len(e) == 1 and (len(inten) == 0 or len(forbidden) == 0):    
        print('File: ', daughter, file=errfile)
        print('No Documented Branches for Parent or Isomer', file=errfile)
    elif len(e) != (len(inten) or len(forbidden)):    
        print('File: ', daughter ,file=errfile)
        print('Length mismatch!', file=errfile)
        print('Length of e: ', len(e),': ',e,file=errfile)
        print('Length of intensity: ', len(inten),': ', inten, file=errfile)
        print('Length of forbidden: ', len(forbidden), ': ', forbidden, file=errfile)

    if i == length:
        print('Successfully Read File', file=errfile)  
    else:
        print('Exception Occured!', file=errfile)
        print('Read terminated early at line', i, file=errfile)

    f.close()
    print('Lines processed', i)

    return length


#Make text file with names of the input files
file = 'cum_py.txt'

file2 = open('betadata.dat', 'w+')
errfile = open('readerr.dat','w+')    

A, element = read_prod_data(file)

for a, elem in zip(A, element):
    isotope_file = a + elem
    
    print(isotope_file)   
    length = data(isotope_file, errfile)
    print('EOF', file=file2)

print('Finished!')
file2.close()
errfile.close()
