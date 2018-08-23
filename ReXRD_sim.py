import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from atomic_form_factors import ff

Avogadro = 6.02214179*10**23
r_electron = 2.8179402894*10**(-13)
evA = 12398.4187526
Henke_path = 'ScatteringFactors\\'

# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
cosd = lambda x: np.cos(x*np.pi/180.)
acosd = lambda x: 180.*np.arccos(x)/np.pi

def parseCIF(filename):
    while True:
        try:
            cif = open(filename, 'r')
            break
        except:
            print "Cannot open the file " + filename
            filename = raw_input("Please re-enter the CIF file and path:  ")
    line = cif.readline()
    lattice = np.zeros(6)
    while not "_symmetry_equiv_pos_as_xyz" in line.lower():
        if "_cell_length_a" in line:
            temp = line.split()
            lattice[0] = float(temp[-1].split('(')[0])
        if "_cell_length_b" in line:
            temp = line.split()
            lattice[1] = float(temp[-1].split('(')[0])
        if "_cell_length_c" in line:
            temp = line.split()
            lattice[2] = float(temp[-1].split('(')[0])
        if "_cell_angle_alpha" in line:
            temp = line.split()
            lattice[3] = float(temp[-1].split('(')[0])
        if "_cell_angle_beta" in line:
            temp = line.split()
            lattice[4] = float(temp[-1].split('(')[0])
        if "_cell_angle_gamma" in line:
            temp = line.split()
            lattice[5] = float(temp[-1].split('(')[0])
        line = cif.readline()
    syms = np.array([np.identity(3), np.identity(3)])
    trans = np.array([np.zeros(3), np.zeros(3)])
    line = cif.readline()
    while True:
        try:
            line = line.translate(None, "'")
            line = line.translate(None, '"')
            temp = line.split(',')
            sym = np.identity(3)
            sym[:,:] = 0.0
            tran = np.zeros(3) 
            lib = {'x' : 0, 'y' : 1, 'z' : 2}
            for i in range(0, 3):
                temp[i] = temp[i].strip()
                temp[i] = temp[i].translate(None, ' ')
                if 'x' in temp[i]:
                    index = temp[i].index('x')
                    if temp[i][index-1] == '-':
                        sym[i,0] = -1.0
                    else:
                        sym[i,0] = 1.0
                if 'y' in temp[i]:
                    index = temp[i].index('y')
                    if temp[i][index-1] == '-':
                        sym[i,1] = -1.0
                    else:
                        sym[i,1] = 1.0
                if 'z' in temp[i]:
                    index = temp[i].index('z')
                    if temp[i][index-1] == '-':
                        sym[i,2] = -1.0
                    else:
                        sym[i,2] = 1.0
                try:
                    num,den = temp[i][-4:].split( '/' )
                    tran[i] = float(num)/float(den)
                except:
                    tran[i] = 0.0
            syms = np.append(syms, [sym], axis=0)
            trans = np.append(trans, [tran], axis=0)
            line = cif.readline()
        except:
            break
    print "Read in " + `len(trans)-2` + " symmetry operators."          
                
    while not "_atom_site_" in line.lower():
        while not "loop_" in line.lower():
            line = cif.readline()
        line = cif.readline()
    count = 0
    while "_atom_site_" in line.lower():
        if "_atom_site_type_symbol" in line.lower():
            elem = count
        if "_atom_site_symmetry_multiplicity" in line.lower():
            mult = count
        if "_atom_site_fract_x" in line.lower():
            x = count
        if "_atom_site_fract_y" in line.lower():
            y = count
        if "_atom_site_fract_z" in line.lower():
            z = count
        if "_atom_site_occupancy" in line.lower():
            occ = count
        count += 1
        line = cif.readline()
    atomlist = []
    while True:
        try:
            temp = line.split()
            s = temp[elem].lower()
            element = ''.join(i for i in s if not i.isdigit() and i not in ['-', '+'])
            multiplicity = int(temp[mult])
            pos = [float(temp[x].split('(')[0]) % 1, float(temp[y].split('(')[0]) % 1, float(temp[z].split('(')[0]) % 1]
            #apply symmetry elements here...
            positions = [pos]
            for i in range(0, len(trans)):
                newpos = np.dot(syms[i], pos) + trans[i]
                newpos = newpos % 1
                count = 0
                for j in range(0, len(positions)):
                    if np.allclose(newpos, positions[j], atol=1e-6):
                        count += 1
                        break
                if count == 0:
                    positions = np.append(positions, [newpos], axis = 0)
            for pos in positions:
                new = Atom()
                new.element = element
                new.xff_from_file()
                new.pos = pos
                new.occ = float(temp[occ].split('(')[0])
                atomlist = np.append(atomlist, new)
            line = cif.readline()
        except:
            break
    cif.close()
    return lattice, atomlist
    
    

def getf(atom):
    while True:
        try:
            atomfile = open(Henke_path + atom + ".nff", 'r')
            break
        except:
            print "Cannot open the atomic scattering factor file for " + atom
            atom = raw_input("Please re-specify the atom name:  ")
    line = atomfile.readline()
    x = []
    f1 = []
    f2 = []
    while line:
        temp = line.split()
        x = np.append(x, float(temp[0]))
        f1 = np.append(f1, float(temp[1]))
        f2 = np.append(f2, float(temp[2]))
        line = atomfile.readline()
    return x, f1, f2

class XFF():
    def __init__(self):
        self.x = np.array([])
        self.f1 = np.array([])
        self.f2 = np.array([])

class Atom(object):
    ''' Atom --> class for storting atom information
    
    Data members:
        element     -- type of the atom
        xff     -- Atomic Scattering factors (energy, f1, f2)
        pos     -- position in fractional coordinates (x, y, z)
        occ     -- occupancy for that site, defaults to 1
    
    Methods:
        xff_from_file   -- Extracts xff data from a .xff file
        xff_interp  -- Returns xff data interpolated over specified energies
    '''                                                                          
    def __init__(self):
        self.element = str()
        self.type = str()
        self.xff = XFF()
        self.pos = [0.0, 0.0, 0.0]
        self.occ = 1.0
    def xff_from_file(self):
        if self.element:
            xff_name = Henke_path + self.element + '.nff'
            #print "obtaining atomic scattering parameters from Henke data for element " + self.element
        else:
            self.element = raw_input("Please specify the atom type (lower case only):  ")
            xff_name = Henke_path + self.element + '.nff'
        try:
            xff_file = open(xff_name, 'r')
            line = xff_file.readline()
        except:
            #print 'Cannot open ' + str(xff_name)
            return False
        self.type = self.element.split('_')[0]
        temp = line.split()
        try:
            float(temp[1])
        except:
            line = xff_file.readline() #throw away header line
        while line:
            temp = line.split()
            self.xff.x = np.append(self.xff.x, float(temp[0]))
            self.xff.f1 = np.append(self.xff.f1, float(temp[1]))
            self.xff.f2 = np.append(self.xff.f2, float(temp[2]))
            line = xff_file.readline()
        xff_file.close()
        return True
    def xff_interp(self, xmesh):
        try:
            indexing = np.argsort(self.xff.x)
            self.xff.x = self.xff.x[indexing]
            self.xff.f1 = self.xff.f1[indexing]
            self.xff.f2 = self.xff.f2[indexing]
            interp_f1 = interpolate.InterpolatedUnivariateSpline(self.xff.x,self.xff.f1, k=1)
            interp_f2 = interpolate.InterpolatedUnivariateSpline(self.xff.x,self.xff.f2, k=1)
            self.xff.f1 = interp_f1(xmesh)
            self.xff.f2 = interp_f2(xmesh)
            self.xff.x = xmesh
        except:
            print 'Cannot interpolate over that range for atom ' + self.name + '.'
    def get_f0(self, Q):
        f0 = ff[self.type.capitalize()][-1]
        for i in xrange(0, 8, 2):
            f0 += ff[self.type.capitalize()][i] *np.exp(-ff[self.type.capitalize()][i+1]*(Q/(4.*np.pi))**2) * np.exp(-0.01 * Q**2 /(4.0*np.pi))
        return f0

def calcS(atomlist, h, k, l, hvs, lattice):
    try:
        lattice = atomlist.lattice
    except:
        pass
    S = np.zeros_like(hvs, type(complex))
    Q = find_Q(lattice, h, k, l)
    for atom in atomlist:
        atom.xff_interp(hvs)
        S[:] += atom.occ*(atom.get_f0(0.0) + (atom.xff.f1[:] - atom.get_f0(0.0)) + 1j*atom.xff.f2[:]) * np.exp(2.0*np.pi*1j*(h*atom.pos[0] + k*atom.pos[1] + l*atom.pos[2]))
    return S

def find_Q(lattice, h, k, l):
    temp1 = 1.0/(1.0 + 2.0*cosd(lattice[3])*cosd(lattice[4])*cosd(lattice[5]) - cosd(lattice[3])**2 - cosd(lattice[4])**2 - cosd(lattice[5])**2)
    A = (float(h)*sind(lattice[3])/lattice[0])**2
    B = (float(k)*sind(lattice[4])/lattice[1])**2
    C = (float(l)*sind(lattice[5])/lattice[2])**2
    D = 2.*float(h)*float(k)*(cosd(lattice[3])*cosd(lattice[4]) - cosd(lattice[5]))/(lattice[0] * lattice[1])
    E = 2.*float(k)*float(l)*(cosd(lattice[4])*cosd(lattice[5]) - cosd(lattice[3]))/(lattice[1] * lattice[2])
    F = 2.*float(l)*float(h)*(cosd(lattice[5])*cosd(lattice[3]) - cosd(lattice[4]))/(lattice[0] * lattice[2])
    return 2.0 * np.pi * np.sqrt(temp1 * (A + B + C + D + E + F))
    
def equiv_peaks(lattice, h, k, l):
    import itertools
    peaklist = []
    Q = find_Q(lattice, h, k, l)
    index_range = np.arange(-8, 9)
    for peak in itertools.product(index_range, repeat=3):
        test_Q = find_Q(lattice, *peak)
        if np.allclose(Q, test_Q):
            peaklist = np.append(peaklist, peak, axis = 0)
    peaklist = np.reshape(peaklist, [-1,3])
    print "Found %i equivlant peaks for %i, %i, %i" %(len(peaklist), h, k, l)
    return peaklist
###############################################################################
''' Example of how to enter atoms by hand if you want to do things that way
cu1 = Atom()
cu1.element = 'cu'
cu1.xff_from_file()
cu1.pos = [0.0, 0.0, 0.0]
'''
if __name__ == "__main__":
    while True:
        name = raw_input("Enter the path and file name of the CIF to use:  ")
        lattice, atomlist = parseCIF(name)
        break
    while True:
        try:
            start_hv = raw_input("Enter low energy (eV):  ")
            start_hv = float(start_hv)
            break
        except:
            print "Not a valid number, try again..."
    while True:
        try:
            stop_hv = raw_input("Enter high energy (eV):  ")
            stop_hv = float(stop_hv)
            break
        except:
            print "Not a valid number, try again..."
    while True:
        try:
            step_hv = raw_input("Enter energy step size (eV):  ")
            step_hv = float(step_hv)
            break
        except:
            print "Not a valid number, try again..."
    hvs = np.arange(start_hv, stop_hv, step_hv)
    peaks = []
    while True:
        try:
            hkl = raw_input("Enter an hkl to use (format h, k, l), blank to end:  ")
            if len(hkl) < 3:
                break
            temp = hkl.split(',')
            h = int(temp[0])
            k = int(temp[1])
            l = int(temp[2])
            peaks.append([h,k,l])
        except:
            print "Could not interpret that h,k,l, try again..."
    i = 0
    for peak in peaks:
        Ycalc = np.zeros_like(hvs)
        h=peak[0]
        k=peak[1]
        l=peak[2]
        equivs = equiv_peaks(lattice, h, k, l)
        for eq in equivs:
            heq = eq[0]
            keq = eq[1]
            leq = eq[2]
            temp = abs(calcS(atomlist, heq, keq, leq, hvs, lattice))**2
            temp2 = Ycalc + temp
            Ycalc = temp2
        plt.figure(i)
        plt.title('(' + str(h) + ', ' + str(k) + ', ' + str(l) + ')')
        plt.plot(hvs, Ycalc, 'r-')
        plt.xlabel('Beamline Energy (eV)')
        plt.ylabel('Peak Intensity (A.U.)')
        i += 1
    plt.show()
