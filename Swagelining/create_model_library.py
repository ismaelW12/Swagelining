from abaqus import *
from abaqusConstants import *
from caeModules import *
import sys
import numpy as np
import csv

###########################################################################
#                                                                         #
#                AUXILLIARY PARAMETERS                                    #
#                                                                         #
###########################################################################

info=sys.stdin.readline()
a=info.split("\\")

#Liner
OD = float(a[0])
WT = float(a[1])
elem_length = float(a[2])
no_el_thru_thkness = int(a[3])
refLength = float(a[4])

#Die
dia_out = float(a[5])
angle = float(a[6])

#Backing pipe
pipeFlag = a[7]
pipeOD = float(a[8])
pipeWT = float(a[9])

#Material
anisotropy = float(a[10])
stressStrainCurve = a[11]
yieldCriterion = a[12]
kFactor = float(a[13])
angle_friction = float(a[14])
ref_stress = float(a[15])
relaxation = float(a[16])

#PE and steel friction
friction = float(a[17])

#Load details
temperature = float(a[18])

#Boundary conditions
BC = a[19]

#File name
filename = a[20]

#Element type
elementType = a[21]

#Auxilliary values
orad = OD / 2
irad = (OD - 2*WT) / 2
pipeID = pipeOD - 2.0 * pipeWT
pipeIrad = 0.5 * pipeID
rad_out = dia_out/2
overlap = (pipeIrad - rad_out) / np.tan(np.radians(angle))
post_rad_out = rad_out + overlap * np.tan(np.radians(angle))
length = (orad - rad_out) / np.tan(np.radians(angle))
rad_in = rad_out + length * np.tan(np.radians(angle))
pre_rad_in = rad_out + (length + overlap) * np.tan(np.radians(angle))
endLength = 100.0 # int(OD / 10.0) * 10.0
L = endLength + refLength + endLength
displ1 = L
displ2 = displ1 + length + overlap + 0.15 * L

###########################################################################
#                                                                         #
#                CREATE PART: LINER                                       #
#                                                                         #
###########################################################################

if elementType == 'axis':

    #Create part
    p1 = mdb.models['Model-1'].Part(name = 'Liner', dimensionality = AXISYMMETRIC, type = DEFORMABLE_BODY)

    #Feature 1: Solid extrude-1
    s1a = mdb.models['Model-1'].ConstrainedSketch(name = '__profile__', sheetSize = 2000.0, gridSpacing = 1.0)
    s1a.rectangle(point1 = (0.0, irad), point2 = (L, orad))
    s1a.ConstructionLine(point1 = (0.0, -100.0), point2 = (0.0, 100.0))
    p1.BaseShell(sketch = s1a)
    f1, e1, d1, v1, c1 = p1.faces, p1.edges, p1.datums, p1.vertices, p1.cells

    #Set 1: set_liner
    p1.Set(faces = f1, name = 'set_liner')

    #Set 2: set_face_pull
    p1.Set(edges = e1.findAt(coordinates = ((0.0, orad - 0.5*WT, 0.0), )), name = 'set_face_pull')

    #Set 3: set_face_fix
    p1.Set(edges = e1.findAt(coordinates = ((L, orad - 0.5*WT, 0.0), )), name = 'set_face_fix')

    #Set 4: set_node_pull
    p1.Set(vertices = v1.findAt(coordinates = ((0.0, orad, 0.0), )), name = 'set_node_pull')

    #Set 5: set_liner_outer
    p1.Set(edges = e1.findAt(coordinates = ((1.0, orad, 0.0), )), name = 'set_liner_outer')

    #Surface 1: surf_liner_outer
    p1.Surface(side1Edges = e1.findAt(coordinates = ((1.0, orad, 0.0), )), name = 'surf_liner_outer')

    #Partition
    s1b = mdb.models['Model-1'].ConstrainedSketch(name = '__profile__', sheetSize = 2000.0, gridSpacing = 1.0)
    s1b.Line(point1 = (endLength, irad), point2 = (endLength, orad))
    s1b.Line(point1 = (endLength + 0.5 * refLength - 0.5 * elem_length, irad), point2 = (endLength + 0.5 * refLength - 0.5 * elem_length, orad))
    s1b.Line(point1 = (endLength + 0.5 * refLength + 0.5 * elem_length, irad), point2 = (endLength + 0.5 * refLength + 0.5 * elem_length, orad))
    s1b.Line(point1 = (endLength + refLength, irad), point2 = (endLength + refLength, orad))
    p1.PartitionFaceBySketch(faces = f1, sketch = s1b)
    f2, e2, d2, v2, c2 = p1.faces, p1.edges, p1.datums, p1.vertices, p1.cells
    p1.Set(faces = f2.findAt(coordinates = ((endLength + 0.5 * refLength, 0.5 * (irad + orad), 0.0), )), name = 'set_output_central')
    p1.Set(faces = f2.findAt(coordinates = ((endLength + 0.5 * refLength, 0.5 * (irad + orad), 0.0), (endLength + 0.5 * refLength - elem_length, 0.5 * (irad + orad), 0.0), (endLength + 0.5 * refLength + elem_length, 0.5 * (irad + orad), 0.0), )), name = 'set_output')

    #Section and assignment
    mdb.models['Model-1'].HomogeneousSolidSection(name = 'section_liner', material = 'mat_liner', thickness = None)
    p1.SectionAssignment(region = p1.sets['set_liner'], sectionName = 'section_liner')

    #Mesh
    #Through wall thickness edge seed
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (0.0, (orad - 0.5 * WT), 0.0)), ), number = no_el_thru_thkness, constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength, (orad - 0.5 * WT), 0.0)), ), number = no_el_thru_thkness, constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength + 0.5 * refLength - 0.5 * elem_length, (orad - 0.5 * WT), 0.0)), ), number = no_el_thru_thkness, constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength + 0.5 * refLength + 0.5 * elem_length, (orad - 0.5 * WT), 0.0)), ), number = no_el_thru_thkness, constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength + refLength, (orad - 0.5 * WT), 0.0)), ), number = no_el_thru_thkness, constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength + refLength + endLength, (orad - 0.5 * WT), 0.0)), ), number = no_el_thru_thkness, constraint = FIXED)

    #Along length edge seed
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (0.5 * endLength, orad, 0.0)), ), number = int(endLength / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength + 0.1 * refLength, orad, 0.0)), ), number = int((0.5 * refLength - 0.5 * elem_length) / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength + 0.5 * refLength, orad, 0.0)), ), number = int(elem_length / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength + 0.9 * refLength, orad, 0.0)), ), number = int((0.5 * refLength - 0.5 * elem_length) / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (endLength + refLength + 0.5 * endLength, orad, 0.0)), ), number = int(endLength / elem_length), constraint = FIXED)

    #Set element type and mesh
    p1.setElementType(elemTypes = (mesh.ElemType(elemCode = CAX4R, elemLibrary = STANDARD),), regions = p1.sets['set_liner'])
    p1.generateMesh()

elif elementType == 'shell':

    #Create part
    p1 = mdb.models['Model-1'].Part(name = 'Liner', dimensionality = THREE_D, type = DEFORMABLE_BODY)

    #Feature 1: Solid extrude-1
    s1a = mdb.models['Model-1'].ConstrainedSketch(name = '__profile__', sheetSize = 2000.0, gridSpacing = 1.0)
    s1a.rectangle(point1 = (irad, -L), point2 = (orad, 0.0))
    s1a.ConstructionLine(point1 = (0.0, -100.0), point2 = (0.0, 100.0))
    p1.BaseSolidRevolve(sketch=s1a, angle=90.0, flipRevolveDirection=OFF)
    f1, e1, d1, v1, c1 = p1.faces, p1.edges, p1.datums, p1.vertices, p1.cells

    #Set 1: set_liner
    p1.Set(cells = c1, name = 'set_liner')

    #Set 2: set_face_pull
    p1.Set(faces = f1.findAt(coordinates = (((orad - 0.5*WT) * np.sin(np.radians(45)), 0.0, (orad - 0.5*WT) * np.cos(np.radians(45))), )), name = 'set_face_pull')

    #Set 3: set_face_fix
    p1.Set(faces = f1.findAt(coordinates = (((orad - 0.5*WT) * np.sin(np.radians(45)), -L, (orad - 0.5*WT) * np.cos(np.radians(45))), )), name = 'set_face_fix')

    #Set 4: set_edge_pull
    p1.Set(edges = e1.findAt(coordinates = ((orad * np.sin(np.radians(45)), 0.0, orad * np.cos(np.radians(45))), )), name = 'set_node_pull')

    #Set 5: set_liner_outer
    p1.Set(faces = f1.findAt(coordinates = ((orad * np.sin(np.radians(45)), -0.5 * L, orad * np.cos(np.radians(45))), )), name = 'set_liner_outer')

    #Set 6: set_face_XSymm
    p1.Set(faces = f1.findAt(coordinates = ((0.0, -0.5 * L, orad - 0.5*WT), )), name = 'set_liner_face_XSymm')

    #Set 7: set_face_ZSymm
    p1.Set(faces = f1.findAt(coordinates = ((orad - 0.5*WT, -0.5 * L, 0.0), )), name = 'set_liner_face_ZSymm')

    #Surface 1: surf_liner_outer
    p1.Surface(side1Faces = f1.findAt(coordinates = ((orad * np.sin(np.radians(45)), -0.5 * L, orad * np.cos(np.radians(45))), )), name = 'surf_liner_outer')

    #Partition
    t1b = p1.MakeSketchTransform(sketchPlane = f1[4], sketchUpEdge = e1[4], sketchPlaneSide = SIDE1, sketchOrientation = RIGHT, origin = (0.0, 0.0, 0.0))
    s1b = mdb.models['Model-1'].ConstrainedSketch(name = '__profile__', sheetSize = 2000.0, gridSpacing = 1.0, transform = t1b)
    p1.projectReferencesOntoSketch(sketch = s1b, filter = COPLANAR_EDGES)
    s1b.Line(point1 = (endLength, irad), point2 = (endLength, orad))
    s1b.Line(point1 = (endLength + 0.5 * refLength - 0.5 * elem_length, irad), point2 = (endLength + 0.5 * refLength - 0.5 * elem_length, orad))
    s1b.Line(point1 = (endLength + 0.5 * refLength + 0.5 * elem_length, irad), point2 = (endLength + 0.5 * refLength + 0.5 * elem_length, orad))
    s1b.Line(point1 = (endLength + refLength, irad), point2 = (endLength + refLength, orad))
    p1.ShellRevolve(sketchPlane=f1[4], sketchUpEdge=e1[4], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s1b, angle=90.0, flipRevolveDirection=ON, keepInternalBoundaries=ON)
    f2, e2, d2, v2, c2 = p1.faces, p1.edges, p1.datums, p1.vertices, p1.cells
    p1.Set(cells = c2.findAt(coordinates = (((orad - 0.5*WT) * np.sin(np.radians(45)), - 0.5 * L, (orad - 0.5*WT) * np.cos(np.radians(45))), )), name = 'set_output_central')
    p1.Set(cells = c2.findAt(coordinates = (((orad - 0.5*WT) * np.sin(np.radians(45)), - 0.5 * L, (orad - 0.5*WT) * np.cos(np.radians(45))), ((orad - 0.5*WT) * np.sin(np.radians(45)), - 0.5 * L + 0.5 * refLength, (orad - 0.5*WT) * np.cos(np.radians(45))), ((orad - 0.5*WT) * np.sin(np.radians(45)), - 0.5 * L - 0.5 * refLength, (orad - 0.5*WT) * np.cos(np.radians(45))))), name = 'set_output')

    #Section and assignment
    mdb.models['Model-1'].HomogeneousShellSection(name = 'section_liner', preIntegrate = OFF, material = 'mat_liner', thicknessType = UNIFORM, thickness = WT, thicknessField = '', nodalThicknessField = '', idealization = NO_IDEALIZATION, poissonDefinition = DEFAULT, thicknessModulus = None, temperature = GRADIENT, useDensity = OFF, integrationRule = SIMPSON, numIntPts = no_el_thru_thkness)
    p1.SectionAssignment(region = p1.sets['set_liner'], sectionName = 'section_liner')

    #Mesh
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (0.0, - 0.5 * endLength, orad), ), ), number = int(endLength / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = ((0.0, -0.5 * L + 0.25 * refLength, orad), (0.0, -0.5 * L - 0.25 * refLength, orad), ))), number = int((0.5 * refLength - 0.5 * elem_length) / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (0.0, -0.5 * L, orad)), ), number = int(elem_length / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (0.0, -L + 0.5 * endLength, orad)), ), number = int(endLength / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (orad * np.sin(np.radians(45)), 0.0, orad * np.cos(np.radians(45)))), ), number = int(orad * np.radians(90) / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (irad * np.sin(np.radians(45)), 0.0, irad * np.cos(np.radians(45)))), ), number = int(orad * np.radians(90) / elem_length), constraint = FIXED)
    p1.seedEdgeByNumber(edges = (e2.findAt(coordinates = (0.0, 0.0, orad - 0.5*WT), ), ), number = 1, constraint = FIXED)

    #Set element type and mesh
    p1.setElementType(elemTypes = (mesh.ElemType(elemCode = SC8R, elemLibrary = STANDARD),), regions = p1.sets['set_liner'])
    p1.generateMesh()

###########################################################################
#                                                                         #
#                CREATE PART: DIE OUTER                                   #
#                                                                         #
###########################################################################

if elementType == 'axis':

    #Create part
    p2 = mdb.models['Model-1'].Part(name = 'Die_Outer', dimensionality = AXISYMMETRIC, type = ANALYTIC_RIGID_SURFACE)

    #Feature 1: Rigid shell-1
    s2 = mdb.models['Model-1'].ConstrainedSketch(name = '__profile__', sheetSize = 2000.0, gridSpacing = 1.0)
    g, v, d = s2.geometry, s2.vertices, s2.dimensions

    s2.Line(point1 = (post_rad_out, length + overlap), point2 = (rad_out, length))
    s2.Line(point1 = (rad_out, length), point2 = (rad_in, 0.0))
    s2.Line(point1 = (rad_in, 0.0), point2 = (pre_rad_in, -overlap))
    s2.FilletByRadius(radius=20.0, curve1=g.findAt((0.5 * (rad_out + rad_in), 0.5 * (length + 0.0))), nearPoint1=(0.5 * (rad_out + rad_in), 0.5 * (length + 0.0)), curve2=g.findAt((0.5 * (post_rad_out + rad_out), 0.5 * (length + overlap +length))), nearPoint2=(0.5 * (post_rad_out + rad_out), 0.5 * (length + overlap +length)))

    s2.ConstructionLine(point1 = (0.0, -100.0), point2 = (0.0, 100.0))
    p2.AnalyticRigidSurf2DPlanar(sketch = s2)
    f2, e2, d2, v2, c2 = p2.faces, p2.edges, p2.datums, p2.vertices, p2.cells

    #Surface 1: surf_die_outer
    p2.Surface(side2Edges = e2, name = 'surf_die_outer')

elif elementType == 'shell':

    #Create part
    p2 = mdb.models['Model-1'].Part(name = 'Die_Outer', dimensionality = THREE_D, type = ANALYTIC_RIGID_SURFACE)

    #Feature 1: Rigid shell-1
    s2 = mdb.models['Model-1'].ConstrainedSketch(name = '__profile__', sheetSize = 2000.0, gridSpacing = 1.0)
    g, v, d = s2.geometry, s2.vertices, s2.dimensions

    s2.Line(point1 = (post_rad_out, length + overlap), point2 = (rad_out, length))
    s2.Line(point1 = (rad_out, length), point2 = (rad_in, 0.0))
    s2.Line(point1 = (rad_in, 0.0), point2 = (pre_rad_in, -overlap))
    s2.FilletByRadius(radius=20.0, curve1=g.findAt((0.5 * (rad_out + rad_in), 0.5 * (length + 0.0))), nearPoint1=(0.5 * (rad_out + rad_in), 0.5 * (length + 0.0)), curve2=g.findAt((0.5 * (post_rad_out + rad_out), 0.5 * (length + overlap +length))), nearPoint2=(0.5 * (post_rad_out + rad_out), 0.5 * (length + overlap +length)))

    s2.ConstructionLine(point1 = (0.0, -100.0), point2 = (0.0, 100.0))
    p2.AnalyticRigidSurfRevolve(sketch = s2)
    f2, e2, d2, v2, c2= p2.faces, p2.edges, p2.datums, p2.vertices, p2.cells

    #Surface 1: surf_die_outer
    p2.Surface(side1Faces = f2, name = 'surf_die_outer')

###########################################################################
#                                                                         #
#             CREATE PART: BACKING PIPE                                   #
#                                                                         #
###########################################################################

if elementType == 'axis':

    #Create part
    p3 = mdb.models['Model-1'].Part(name = 'Backing_Pipe', dimensionality = AXISYMMETRIC, type = ANALYTIC_RIGID_SURFACE)
    
    #Feature 1: Rigid shell-1
    s3 = mdb.models['Model-1'].ConstrainedSketch(name = '__profile__', sheetSize = 2000.0, gridSpacing = 1.0)
    s3.Line(point1 = (pipeIrad, length + overlap + 1.2 * L), point2 = (pipeIrad, length + elem_length + overlap - 0.25 * overlap))
    s3.ConstructionLine(point1 = (0.0, -100.0), point2 = (0.0, 100.0))
    p3.AnalyticRigidSurf2DPlanar(sketch = s3)
    f3, e3, d3, v3, c3 = p3.faces, p3.edges, p3.datums, p3.vertices, p3.cells
    
    #Surface 1: surf_backing_pipe
    p3.Surface(side2Edges = e3, name = 'surf_backing_pipe')

elif elementType == 'shell':

    #Create part
    p3 = mdb.models['Model-1'].Part(name = 'Backing_Pipe', dimensionality = THREE_D, type = ANALYTIC_RIGID_SURFACE)
    
    #Feature 1: Rigid shell-1
    s3 = mdb.models['Model-1'].ConstrainedSketch(name = '__profile__', sheetSize = 2000.0, gridSpacing = 1.0)
    s3.Line(point1 = (pipeIrad, length + overlap + 1.2 * L), point2 = (pipeIrad, length + elem_length + overlap - 0.25 * overlap))
    s3.ConstructionLine(point1 = (0.0, -100.0), point2 = (0.0, 100.0))
    p3.AnalyticRigidSurfRevolve(sketch = s3)
    f3, e3, d3, v3, c3 = p3.faces, p3.edges, p3.datums, p3.vertices, p3.cells

    #Surface 1: surf_backing_pipe
    p3.Surface(side2Faces = f3, name = 'surf_backing_pipe')

###########################################################################
#                                                                         #
#                CREATE MATERIAL                                          #
#                                                                         #
###########################################################################

#Create material
m1 = mdb.models['Model-1'].Material(name = 'mat_liner')

#Create empty arrays for material properties
temperatures       = np.empty(0)
elastic_moduli     = np.empty(0)
poissons           = np.empty(0)
angles_friction    = np.empty(0)
angles_dilation    = np.empty(0)
ref_stresses       = np.empty(0)
stresses           = np.empty(0)
strains_0          = np.empty(0)
strains_23         = np.empty(0)
strains_60         = np.empty(0)

#Read stress-strain data from csv
with open('stress-strain.csv') as file:
    reader = csv.DictReader(file)
    for row in reader:
        stresses   = np.append(  stresses, float(row['Stress']))
        strains_0  = np.append( strains_0, float(row['Strain_0']))
        strains_23 = np.append(strains_23, float(row['Strain_23']))
        strains_60 = np.append(strains_60, float(row['Strain_60']))

#Read material data from csv
with open('material_data.csv') as file:
    reader = csv.DictReader(file)
    for row in reader:
        temperatures    = np.append(   temperatures, float(row['temperature']))
        elastic_moduli  = np.append( elastic_moduli, float(row['elastic_modulus']))
        poissons        = np.append(       poissons, float(row['poisson']))
        angles_friction = np.append(angles_friction, float(row['angle_friction']))
        angles_dilation = np.append(angles_dilation, float(row['angle_dilation']))
        ref_stresses    = np.append(   ref_stresses, float(row['reference_stress']))

#Interpolate for material properties based on temperature
elastic_modulus = np.interp(temperature, temperatures, elastic_moduli)
poisson         = np.interp(temperature, temperatures, poissons)
#angle_friction  = np.interp(temperature, temperatures, angles_friction)
angle_dilation  = np.interp(temperature, temperatures, angles_dilation)
#ref_stress      = np.interp(temperature, temperatures, ref_stresses)

#Calculate drucker prager hardening
trueStress = stresses
trueStrain = np.zeros(trueStress.size)
for i in range(strains_0.size):
    trueStrain[i] = np.interp(temperature, [0.0, 23.0, 60.0], [strains_0[i], strains_23[i], strains_60[i]])
truePlasticStrain = trueStrain - ref_stress / elastic_modulus

stressStrainCurveTable = []
if stressStrainCurve == 'EPP':
    abqStrain = np.array([0.0, 1.0])
    abqStress = np.array([ref_stress, 1.01 * ref_stress])
    stressStrainCurveTable.append((abqStress[0], abqStrain[0]))
    stressStrainCurveTable.append((abqStress[1], abqStrain[1]))
elif stressStrainCurve == 'EP':
    abqStrain=np.concatenate((np.arange(0.0, 0.1, 0.005),np.arange(0.1, 1.05, 0.05)))
    abqStress=np.interp(abqStrain, truePlasticStrain, trueStress)
    for i in range(abqStress.size):
        stressStrainCurveTable.append((abqStress[i], abqStrain[i]))

#Set material properties
if anisotropy == 1.0:
    #m1.Elastic(moduli=INSTANTANEOUS, table=((elastic_modulus, poisson), ))
    m1.Elastic(moduli=INSTANTANEOUS, table=((1200.0, 0.45), ))
else:
    m1.Elastic(type=ENGINEERING_CONSTANTS, table=((anisotropy * elastic_modulus, elastic_modulus, anisotropy * elastic_modulus, poisson, anisotropy * poisson, poisson, 0.5 * elastic_modulus / (1 + poisson), 0.5 * elastic_modulus / (1 + poisson), 0.5 * elastic_modulus / (1 + poisson)), ))
if relaxation == 1.0:
    m1.Viscoelastic(domain=TIME, time=PRONY, table=((0.999, 0.999, 7.0), ))
elif relaxation == 2.0:
    m1.Viscous(table=((0.001, 2.4, -0.5, 0.80), ))
if stressStrainCurve != 'E' and yieldCriterion == 'M':
    #None
    #m1.Plastic(table = stressStrainCurveTable)
    m1.Plastic(table = ((23.0, 0.0), (27.979253112, 1.0), ))
    if relaxation == 1.0:
        None
        #m1.Viscous(law=TIME, table=((1.0, 1.0, -0.5, 0.5), ))
elif stressStrainCurve != 'E' and yieldCriterion == 'DP':
    m1.DruckerPrager(table = ((angle_friction, kFactor, angle_dilation), ))
    m1.druckerPrager.DruckerPragerHardening(table = stressStrainCurveTable, type = TENSION)
    #m1.druckerPrager.DruckerPragerCreep(table = ((2.0E-04, 1.5, -0.1), ))

###########################################################################
#                                                                         #
#                CREATE ASSEMBLY                                          #
#                                                                         #
###########################################################################

#Create assembly
a = mdb.models['Model-1'].rootAssembly
a.deleteAllFeatures()
a.DatumCsysByThreePoints(coordSysType = CYLINDRICAL, origin = (0.0, 0.0, 0.0), point1 = (1.0, 0.0, 0.0), point2 = (0.0, 0.0, -1.0))
a_liner = a.Instance(name = 'Liner-1', part = p1, dependent = ON)
a_die_outer = a.Instance(name = 'Die_Outer-1', part = p2, dependent = ON)
if pipeFlag == 'True':
    a_backing_pipe = a.Instance(name = 'Backing_Pipe-1', part = p3, dependent = ON)

#Rotate and translate instances
if elementType == 'axis':
    a.rotate(instanceList = ('Liner-1', ), angle = 270.0, axisPoint = (0.0, 0.0, 0.0), axisDirection = (0.0, 0.0, 1.0))

#Reference point
RP1 = a.ReferencePoint(point = (0.0, 0.0, 0.0))
a.Set(name = 'geom_RP1', referencePoints = (a.referencePoints[RP1.id], ))
if pipeFlag == 'True':
    RP2 = a.ReferencePoint(point = (0.0, 0.0, 0.0))
    a.Set(name = 'geom_RP2', referencePoints = (a.referencePoints[RP2.id], ))
#
#Orientation
mdb.models['Model-1'].parts['Liner'].MaterialOrientation(region = a_liner.sets['set_liner'], orientationType = GLOBAL, axis = AXIS_3, additionalRotationType = ROTATION_NONE, localCsys = None, fieldName = '', stackDirection = STACK_3)

#Create steps
if relaxation == 1.0:
    mdb.models['Model-1'].ViscoStep(name = 'Step-1', description = 'Insertion', previous = 'Initial', nlgeom = ON, timePeriod = 100.0, initialInc = 0.1, maxInc = 1.0, maxNumInc = 1000, minInc = 1.0E-06, cetol = 0.0001, stabilizationMagnitude=0.0002 ,stabilizationMethod=DISSIPATED_ENERGY_FRACTION, continueDampingFactors=False, adaptiveDampingRatio=0.05)
    mdb.models['Model-1'].ViscoStep(name = 'Step-2', description = 'Reversion', previous = 'Step-1', nlgeom = ON, timePeriod = 5000.0, initialInc = 0.1, maxInc = 25.0, maxNumInc = 1000, minInc = 1.0E-06, cetol = 0.0001, stabilizationMagnitude=0.0002 ,stabilizationMethod=DISSIPATED_ENERGY_FRACTION, continueDampingFactors=False, adaptiveDampingRatio=0.05)    
elif relaxation == 2.0:
    mdb.models['Model-1'].ViscoStep(name = 'Step-1', description = 'Insertion', previous = 'Initial', nlgeom = ON, timePeriod = 0.01, initialInc = 0.0001, maxInc = 0.0001, maxNumInc = 5000, minInc = 1.0E-06, cetol = 0.01, stabilizationMagnitude=1.0E-07 ,stabilizationMethod=DISSIPATED_ENERGY_FRACTION, continueDampingFactors=False, adaptiveDampingRatio=0.05)
    mdb.models['Model-1'].ViscoStep(name = 'Step-2', description = 'Reversion', previous = 'Step-1', nlgeom = ON, timePeriod = 3600.0, initialInc = 0.1, maxInc = 25.0, maxNumInc = 5000, minInc = 1.0E-06, cetol = 0.01, stabilizationMagnitude=1.0E-07 ,stabilizationMethod=DISSIPATED_ENERGY_FRACTION, continueDampingFactors=False, adaptiveDampingRatio=0.05)    
else:
    mdb.models['Model-1'].StaticStep(name = 'Step-1', description = 'Insertion', previous = 'Initial', nlgeom = ON, timePeriod = 5.0, initialInc = 0.01, maxInc = 0.01, maxNumInc = 500, minInc = 1.0E-06)
    mdb.models['Model-1'].StaticStep(name = 'Step-2', description = 'Reversion', previous = 'Step-1', nlgeom = ON, timePeriod = 5.0, initialInc = 0.01, maxInc = 0.01, maxNumInc = 500, minInc = 1.0E-06)

#Interaction property
int_prop1 = mdb.models['Model-1'].ContactProperty('IntProp-1')
int_prop1.TangentialBehavior(formulation = PENALTY, table = ((friction, ), ), fraction = 0.005)
int_prop1.NormalBehavior(allowSeparation = ON, constraintEnforcementMethod = DEFAULT, pressureOverclosure = HARD)
int_prop2 = mdb.models['Model-1'].ContactProperty('IntProp-2')
int_prop1.TangentialBehavior(formulation = PENALTY, table = ((0.0, ), ), fraction = 0.005)
int_prop2.NormalBehavior(allowSeparation = ON, constraintEnforcementMethod = DEFAULT, pressureOverclosure = HARD)

#Interaction
mdb.models['Model-1'].SurfaceToSurfaceContactStd(name = 'Int-1', createStepName = 'Initial', master = a_die_outer.surfaces['surf_die_outer'], slave = a_liner.surfaces['surf_liner_outer'],  interactionProperty='IntProp-1',  sliding=FINITE, enforcement=NODE_TO_SURFACE, thickness=OFF, surfaceSmoothing=NONE, adjustMethod=NONE, smooth=0.2, initialClearance=OMIT, datumAxis=None, clearanceRegion=None)   
if pipeFlag == 'True':
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name = 'Int-2', createStepName = 'Initial', master = a_backing_pipe.surfaces['surf_backing_pipe'], slave = a_liner.surfaces['surf_liner_outer'],  interactionProperty='IntProp-1',  sliding=FINITE, enforcement=NODE_TO_SURFACE, thickness=OFF, surfaceSmoothing=NONE, adjustMethod=NONE, smooth=0.2, initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

#Constraint
mdb.models['Model-1'].RigidBody(name = 'Constraint-1', refPointRegion = a.sets['geom_RP1'], surfaceRegion = a_die_outer.surfaces['surf_die_outer'])
if pipeFlag == 'True':
    mdb.models['Model-1'].RigidBody(name = 'Constraint-2', refPointRegion = a.sets['geom_RP2'], surfaceRegion = a_backing_pipe.surfaces['surf_backing_pipe'])

#Boundary conditions
mdb.models['Model-1'].EncastreBC(name = 'BC-1', createStepName = 'Initial', region = a.sets['geom_RP1'])
if pipeFlag == 'True':
    mdb.models['Model-1'].EncastreBC(name = 'BC-2', createStepName = 'Initial', region = a.sets['geom_RP2'])
mdb.models['Model-1'].DisplacementBC(name = 'BC-3', createStepName = 'Initial', region = a_liner.sets['set_node_pull'], u2 = SET)
if elementType == 'shell':
    mdb.models['Model-1'].DisplacementBC(name = 'BC-XSymm', createStepName = 'Initial', region = a_liner.sets['set_liner_face_XSymm'], u1 = SET)
    mdb.models['Model-1'].DisplacementBC(name = 'BC-ZSymm', createStepName = 'Initial', region = a_liner.sets['set_liner_face_ZSymm'], u3 = SET)
mdb.models['Model-1'].boundaryConditions['BC-3'].setValuesInStep(stepName = 'Step-1', u2 = displ1)
if BC == 'base':
    mdb.models['Model-1'].boundaryConditions['BC-3'].setValuesInStep(stepName = 'Step-2', u2 = displ2)
elif BC == 'alt':
    mdb.models['Model-1'].boundaryConditions['BC-3'].deactivate('Step-2')
    mdb.models['Model-1'].DisplacementBC(name='BC-4', createStepName='Step-2', region=a_liner.sets['set_face_fix'], u1=UNSET, u2=SET, ur3=UNSET, amplitude=UNSET, fixed=ON, distributionType=UNIFORM, fieldName='', localCsys=None)

#Output requests
sectionPointList = []
for i in range(no_el_thru_thkness):
    sectionPointList.append(i+1)
mdb.models['Model-1'].FieldOutputRequest(name = 'F-Output-1', createStepName = 'Step-1', timeInterval = 0.01, timeMarks = OFF, variables = ('COORD', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', ))
mdb.models['Model-1'].FieldOutputRequest(name = 'F-Output-2', createStepName = 'Step-1', timeInterval = 0.01, timeMarks = OFF, variables = ('S', 'LE', 'PE', 'PEEQ', 'PEMAG', ), sectionPoints = sectionPointList)

###########################################################################
#                                                                         #
#                CREATE INPUT FILE AND SAVE MODEL                         #
#                                                                         #
###########################################################################

#Job
job = mdb.Job(name = filename, model = 'Model-1')
job.writeInput()

#Save model
mdb.saveAs(pathName = filename)