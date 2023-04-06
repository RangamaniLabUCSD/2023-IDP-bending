import pymem3dg as dg
import pymem3dg.util as util
import pymem3dg.broilerplate as dg_broil
import numpy as np
from functools import partial
import pdb


####################################################
#                 Initialize paths                #
####################################################
outputDir = "./"

####################################################
#            Initialize input geometry             #
####################################################
""" Built-in construction """
Face, Vertex = dg.getIcosphere(1, 5)
# Vertex = dgu.sphericalHarmonicsPerturbation(Vertex, 5, 6, 0.1)

####################################################
#                 Parameters                       #
####################################################
p = dg.Parameters()

# p.temperature = 0
p.boundary.shapeBoundaryCondition = "none"


p.protein.proteinInteriorPenalty = (
    0  # Constraint to keep protein concentration between 0 and 1
)
p.boundary.proteinBoundaryCondition = "none"
p.variation.isProteinVariation = True

p.proteinMobility = 0.1  # μm^2/s?
# p.adsorption.epsilon = -1 * Kb / R_bar**2
p.adsorption.epsilon = -1e-2
# p.aggregation.chi = 20 * Kb / R_bar**2
# p.aggregation.chi = -0.06
# p.dirichlet.eta = 0.1

p.variation.isProteinConservation = False

p.variation.isShapeVariation = True
p.variation.geodesicMask = -1  # No geodesic mask

# Deviatoric modulus...
p.bending.Kd = 0
p.bending.Kdc = 0

# Bending modulus...
Kb = 8.22e-5  # Bending modulus μm*nN corresponds to ~20 KT
p.bending.Kb = Kb
p.bending.Kbc = 0  # 8.22e-4 #DEFINITION OF LARGE AND SMALL VALUE
# Protein spontaneous curvature
p.bending.H0c = -12


# p.tension.form = partial(dg_broil.constantSurfaceTensionModel, tension=1)  # nN/μΜ

p.tension.form = partial(
    dg_broil.preferredAreaSurfaceTensionModel,
    modulus=0.5,
    preferredArea=4 * np.pi * np.power(1, 2),
)

p.osmotic.form = partial(
    dg_broil.preferredVolumeOsmoticPressureModel,
    preferredVolume=0.8 * (4 / 3 * np.pi * np.power(1, 3)),
    reservoirVolume=0,
    strength=0.01,
)

p.selfAvoidance.d = 1e-5
p.selfAvoidance.mu = 1e-7
p.selfAvoidance.p = 0       # Update every step
p.selfAvoidance.n = 2


p.external.form = partial(
    dg_broil.prescribeGaussianPointForce, Kf=-40* Kb, std=0.01, tau=8
)


if __name__ == "__main__":
    ####################################################
    #                 System                           #
    ####################################################
    """System construction"""
    # Set vertex 0 to be the notable vertex
    notableVertex = np.full(np.shape(Vertex)[0], False)
    notableVertex[0] = True
    proteinDensity = np.full(np.shape(Vertex)[0], 0.1)
    velocity = np.zeros(np.shape(Vertex))

    geometry = dg.Geometry(
        faceMatrix=Face,
        vertexMatrix=Vertex,
        referenceVertexMatrix=Vertex,
        notableVertex=notableVertex,
    )

    g = dg.System(
        geometry=geometry,
        parameters=p,
        velocity=velocity,
        proteinDensity=proteinDensity,
    )

    ####################################################
    #                 Mesh processor                   #
    ####################################################
    g.meshProcessor.meshMutator.mutateMeshPeriod = 10
    g.meshProcessor.meshMutator.isShiftVertex = True
    g.meshProcessor.meshMutator.flipNonDelaunay = True
    # g.meshProcessor.meshMutator.splitLarge = True
    g.meshProcessor.meshMutator.splitFat = True
    g.meshProcessor.meshMutator.splitSkinnyDelaunay = True
    g.meshProcessor.meshMutator.splitCurved = True
    g.meshProcessor.meshMutator.minimumEdgeLength = 0.0001
    g.meshProcessor.meshMutator.maximumEdgeLength = 0.2
    g.meshProcessor.meshMutator.curvTol = 0.005
    g.meshProcessor.meshMutator.collapseSkinny = False 
    g.meshProcessor.meshMutator.collapseSmall = False
    g.meshProcessor.meshMutator.collapseFlat = True 
    g.meshProcessor.meshMutator.targetFaceArea = 0.00005
    g.meshProcessor.meshMutator.isSmoothenMesh = True

    g.initialize(ifMutateMesh=False, ifMute=False)

    # pdb.set_trace()

    ###################################################
    #          Time integration / Optimization
    ####################################################
    """ Integrator setups (essential) """
    h = 0.25
    T = 50000
    eps = 1e-8
    tSave = 100 * h

    # pymem3dg._core.Euler(system: mem3dg::solver::System, characteristicTimeStep: float, totalTime: float, savePeriod: float, tolerance: float, outputDirectory: str, frame: int = 0)
    """ Integrator construction """
    fe = dg.Euler(
        system=g,
        characteristicTimeStep=h,
        totalTime=T,
        savePeriod=tSave,
        tolerance=eps,
        outputDirectory=outputDir,
        frame=0,
    )

    """ Integrator setups (optional) """
    fe.isBacktrack = False
    # fe.ifAdaptiveStep = True
    fe.ifPrintToConsole = True
    fe.ifOutputTrajFile = True
    # fe.ifOutputMeshFile = True
    fe.integrate()
