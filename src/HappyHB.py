class Atom:

    def __init__(self, AtomPdbNumber, AtomTypeExtended, ResidueName, Chain, ResidueNumber, Xcoordinate, Ycoordinate, Zcoordinate, AtomType):
        self.AtomPdbNumber = AtomPdbNumber
        self.AtomTypeExtended = AtomTypeExtended
        self.ResidueName = ResidueName
        self.Chain = Chain
        self.ResidueNumber = ResidueNumber
        self.Xcoordinate = Xcoordinate
        self.Ycoordinate = Ycoordinate
        self.Zcoordinate = Zcoordinate
        self.AtomType = AtomType

        self.Hybridisation = None
    
        self.AtomBonds = []
        self.BondPartners = []

    def GetAtomPdbNumber(self):
        return self.AtomPdbNumber

    def GetAtomTypeExtended(self):
        return self.AtomTypeExtended

    def GetResidueName(self):
        return self.ResidueName

    def GetChain(self):
        return self.Chain

    def GetResidueNumber(self):
        return self.ResidueNumber

    def GetXcoordinate(self):
        return self.Xcoordinate

    def GetYcoordinate(self):
        return self.Ycoordinate

    def GetZcoordinate(self):
        return self.Zcoordinate

    def GetAtomType(self):
        return self.AtomType

    def	GetVdWRadius(self):
        if self.AtomType == 'C':
            return float(1.70)
        elif self.AtomType == 'O':
            return float(1.52)
        elif self.AtomType == 'N':
            return float(1.55)
        elif self.AtomType == 'H':
            return float(1.09)
        elif self.AtomType == 'S':
            return float(1.80)
        elif self.AtomType == 'P':
            return float(1.80)
        elif self.AtomType == 'F':
            return float(1.47)
        elif self.AtomType == 'Br':
            return float(1.85)
        elif self.AtomType == 'I':
            return float(1.98)
        elif self.AtomType == 'X':
            return float(2.00)
        else:    #otherwise - if the above conditions don't satisfy(are not True)
            print("Atom's VdV radius not found")
            return False

    def AddBondPartner(self, Partner):
        self.BondPartners.append(Partner)

    def GetBondPartner(self):
        return self.BondPartners

    def	SetHybridisation(self,Hybridisation):
        self.Hybridisation = Hybridisation

    def	GetHybridisation(self):
        return self.Hybridisation

    def IsDonor(self):
        DonorAtoms = ['H']
       
        if self.AtomType in DonorAtoms:
            for PartnerinHB in self.BondPartners:
                
                if PartnerinHB.AtomType == 'N' or PartnerinHB.AtomType == 'O' or PartnerinHB.AtomType == 'S':
                    return True
                else:
                    return False 
        else:
            return False

    def IsAcceptor(self):
        AcceptorAtoms = ['O','N','S']
        
        if self.AtomType in AcceptorAtoms:
            if self.Hybridisation == 1 or self.Hybridisation == 2:
                return True
            else:
                return False

        else:
            return False

class Residue():
    AtomsInResidue =[]

    def __init__(self, AtomList):
        self.AtomsInResidue.append(AtomList)

def ReadPbdFile(PdbFile):
    AtomsList = []

    with open(PdbFile) as f:
        Lines = f.readlines()

    for Line in Lines:
        LineType = Line[0:6].replace(" ", "")

        if LineType == 'ATOM' or LineType == 'HETATM':
            AtomPdbNumber = int(Line[6:11].replace(" ", ""))
            AtomTypeExtended = Line[12:16].replace(" ", "")
            ResidueName = Line[17:20].replace(" ", "")
            Chain = Line[21].replace(" ", "")
            ResidueNumber = int(Line[22:31].replace(" ", ""))
            Xcoordinate = float(Line[30:38].replace(" ", ""))
            Ycoordinate = float(Line[38:46].replace(" ", ""))
            Zcoordinate = float(Line[46:54].replace(" ", ""))
            AtomType = Line[76:78].replace(" ", "")

            AtomsList.append(Atom(AtomPdbNumber, AtomTypeExtended, ResidueName, Chain, ResidueNumber, Xcoordinate, Ycoordinate, Zcoordinate, AtomType))

    return AtomsList

def AtomsDistance(Atom1, Atom2):
    Distance = ((Atom1.Xcoordinate - Atom2.Xcoordinate)**2 + (Atom1.Ycoordinate - Atom2.Ycoordinate)**2 + (Atom1.Zcoordinate - Atom2.Zcoordinate)**2) ** (1./2.)

    return Distance

def Bonded(Atom1, Atom2):

    if (AtomsDistance(Atom1, Atom2) < (Atom1.GetVdWRadius() + Atom2.GetVdWRadius()) * 0.528):
        return True
    else:
        return False    

def CreateConnectivityMatrix(AtomList):
    i = 0
    for Atom1Instancece in AtomList:
        i = i + 1
        # print(i)
        for Atom2Instancece in AtomList[i:]:
            
            
            if Bonded(Atom1Instancece, Atom2Instancece) == True:
                Atom1Instancece.AddBondPartner(Atom2Instancece)
                Atom2Instancece.AddBondPartner(Atom1Instancece)

def AssignHybridisation(AtomList):

    for Atom1Instancece in AtomList:
        
        NumberOfPartners = len(Atom1Instancece.GetBondPartner())
        Atom1Instancece.SetHybridisation(NumberOfPartners)

# def ClosestNeighbour(Atom1, AtomList):
#     AtomNeighbours={}

#     for Neighbour in AtomList:
#         Distance = AtomsDistance(Atom1, Neighbour)
#         AtomNeighbours[Neighbour] = Distance

#     ClosestSorted = sorted(AtomNeighbours.items(), key=lambda x: x[1], reverse=False)

#     for ClosestFound, Distance in ClosestSorted:

#         if Atom1.GetResidueNumber() != ClosestFound.GetResidueNumber():
#             return ClosestFound

#     return False

def InvolvedInHB(Atom1, AtomList):

    Closest = None
    AtomNeighbours={}

    for Neighbour in AtomList:
        Distance = AtomsDistance(Atom1, Neighbour)
        AtomNeighbours[Neighbour] = Distance

    ClosestSorted = sorted(AtomNeighbours.items(), key=lambda x: x[1], reverse=False)

    for ClosestFound, Distance in ClosestSorted:

        if Atom1.GetResidueNumber() != ClosestFound.GetResidueNumber():
            Closest = ClosestFound
            break

    if Closest == None:
        return False

    if Atom1.IsDonor():
        if Closest.IsAcceptor() and (AtomsDistance(Atom1, Closest) < float(3.00)):
            return True
        else:
            return False
    if Atom1.IsAcceptor():
        if Closest.IsDonor() and (AtomsDistance(Atom1, Closest) < float(3.00)):
            return True
        else:
            return False
    else:
        return False

def GenerateResidueList(AtomList):
    ResidueList = []
    tmpResidueNumber = AtomList[0].GetResidueNumber()
    for Atom1 in AtomList:
        tmpResidueList = []
        if tmpResidueNumber == Atom1.GetResidueNumber():
            tmpResidueList.append(Atom1)
        else:
            ResidueList.append(Residue(tmpResidueList))
            tmpResidueNumber = Atom1.GetResidueNumber()

    return ResidueList

def AddWater(Atom1, AtomList):
    partner = Atom1.GetBondPartner()[0]
    vector = [Atom1.GetXcoordinate() - partner.GetXcoordinate(), Atom1.GetYcoordinate() - partner.GetYcoordinate(), Atom1.GetZcoordinate() - partner.GetZcoordinate()]
    
    Xcoordinate = Atom1.GetXcoordinate() + 2 * vector[0]
    Ycoordinate = Atom1.GetYcoordinate() + 2 * vector[1]
    Zcoordinate = Atom1.GetZcoordinate() + 2 * vector[2]

    print(round(Xcoordinate, 3), round(Ycoordinate, 3), round(Zcoordinate, 3))

    AtomList.append(Atom('1', 'X', 'X', 'X', 'X', Xcoordinate, Ycoordinate, Zcoordinate, 'X'))

    return AtomList

def main(file):

    AtomList = ReadPbdFile(file)

    ResidueList = GenerateResidueList(AtomList)

    CreateConnectivityMatrix(AtomList)
    AssignHybridisation(AtomList)

    for Atom1Instancece in AtomList:
        if Atom1Instancece.IsDonor() and not InvolvedInHB(Atom1Instancece,AtomList) and Atom1Instancece.GetResidueName() == 'UNK':
            AtomList = AddWater(Atom1Instancece, AtomList)

    CreateConnectivityMatrix(AtomList)
    AssignHybridisation(AtomList)
            
    n = 0
    for AtomInstance1 in AtomList:
        if AtomInstance1.GetAtomType() == 'X':
            if AtomInstance1.GetHybridisation() != 0:
                n = n +1
    print('Found', n, 'Unhappy HB donors')
    return n
