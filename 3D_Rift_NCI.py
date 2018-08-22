
# coding: utf-8

# # 3D Lithospheric Model

# In[1]:


import UWGeodynamics as GEO
import os

# get resolution factor for NCI testing
try:
    factor = int(os.environ["UW_RESFACTOR"])
except:
    factor = 1

base_resolution = [128,64,64]
resolution      = [i * factor for i in base_resolution]


# In[2]:


u = GEO.UnitRegistry


# In[3]:


# Characteristic values of the system
half_rate = 1.8 * u.centimeter / u.year
model_length = 500e3 * u.meter
model_width = 500e3 * u.meter
surfaceTemp = 273.15 * u.degK
baseModelTemp = 1603.15 * u.degK
bodyforce = 3370 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT


# In[4]:


Model = GEO.Model(elementRes=resolution, 
                  minCoord=(0. * u.kilometer, 0. * u.kilometer, -160. * u.kilometer), 
                  maxCoord=(500. * u.kilometer, 500. * u.kilometer, 20. * u.kilometer),
                  periodic=[False, True, False],
                  gravity=(0.0, 0.0, -9.81 * u.meter / u.second**2))


# In[5]:


Model.outputDir="outputs"


Model.maxViscosity = 5e23 * u.pascal * u.second
Model.minViscosity = 1e19 * u.pascal * u.second
Model.stressLimiter = 300. * u.megapascal

# In[6]:


Model.diffusivity = 1e-6 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)


# In[7]:


air        = Model.add_material(name="Air", shape=GEO.shapes.Layer(top=Model.top, bottom=0 * u.kilometer))
crust = Model.add_material(name="Crust", shape=GEO.shapes.Layer(top=air.bottom, bottom=-40 * u.kilometer))
mantleLithosphere = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer(top=crust.bottom, bottom=-100 * u.kilometer))
mantle     = Model.add_material(name="Mantle", shape=GEO.shapes.Layer(top=mantleLithosphere.bottom, bottom=Model.bottom))
sediment   = Model.add_material(name="Sediment")

# In[8]:


air.capacity = 100. * u.joule / (u.kelvin * u.kilogram)


# In[9]:


air.density         = 1. * u.kilogram / u.metre**3
crust.density  = GEO.LinearDensity(reference_density=2800. * u.kilogram / u.metre**3)
mantleLithosphere.density  = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3)
mantle.density      = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3)
sediment.density = 2300. * u.kilogram / u.metre**3


# In[10]:


crust.radiogenicHeatProd = 0.7 * u.microwatt / u.meter**3
sediment.radiogenicHeatProd = 0.6 * u.microwatt / u.meter**3


# In[11]:


rh = GEO.ViscousCreepRegistry()


# In[12]:


air.viscosity                = 1e19 * u.pascal * u.second
crust.viscosity         = 1.0 * rh.Gleason_and_Tullis_1995
mantleLithosphere.viscosity  = 5.0 * rh.Karato_and_Wu_1990
mantle.viscosity             = rh.Karato_and_Wu_1990
sediment.viscosity           = rh.Gleason_and_Tullis_1995


# In[13]:


pl = GEO.PlasticityRegistry()


# In[14]:


crust.plasticity         = pl.Huismans_et_al_2011_Crust
mantleLithosphere.plasticity  = pl.Huismans_et_al_2011_Crust
mantle.plasticity             = pl.Huismans_et_al_2011_Crust
sediment.plasticity             = pl.Huismans_et_al_2011_Crust
crust.plasticity.epsilon1 = 0.01
mantleLithosphere.plasticity.epsilon1 = 0.01
mantle.plasticity.epsilon1 = 0.01
sediment.plasticity.epsilon1 = 0.01
crust.plasticity.epsilon2 = 1.0
mantleLithosphere.plasticity.epsilon2 = 1.0
mantle.plasticity.epsilon2 = 1.0
sediment.plasticity.epsilon2 = 1.0

# ## Passive Tracers

# In[15]:


import numpy as np

xp = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), 100)
yp = np.linspace(GEO.nd(Model.minCoord[1]), GEO.nd(Model.maxCoord[1]), 100)

xp, yp = np.meshgrid(xp, yp)
xp = xp.flatten()
yp = yp.flatten()
zp = np.zeros(xp.shape)

surface_tracers = Model.add_passive_tracers(name="Surface", vertices=[xp, yp, zp])
moho_tracers = Model.add_passive_tracers(name="Moho", vertices=[xp, yp, zp+GEO.nd(mantleLithosphere.top)])


# ## Temperature Boundary Condition

# In[16]:


Model.set_temperatureBCs(top=293.15 * u.degK, 
                         bottom=1603.15 * u.degK, 
                         indexSets=[(mantle.indices, 1603.15 * u.degK),
                                    (air.indices, 293.15 * u.degK )])



# In[17]:

velFunc = -Model.y / GEO.nd(Model.maxCoord[1]) * GEO.nd(2.0 * u.centimeter / u.year) + GEO.nd(2.5 * u.centimetre / u.year)

Model.set_velocityBCs(left=[-2.5 * u.centimetre / u.year, None, None],
                      right=[2.5 * u.centimetre / u.year, None, None],
                      bottom=GEO.LecodeIsostasy(reference_mat=mantle.index,
                                                average=False))


# In[18]:


import numpy as np

def gaussian(xx, centre, width):
    return ( np.exp( -(xx - centre)**2 / width ))

maxDamage = 0.25
Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], (GEO.nd(Model.maxCoord[0] - Model.minCoord[0])) / 2.0, GEO.nd(5.0 * u.kilometer))
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,2], GEO.nd(-35. * u.kilometer) , GEO.nd(5.0 * u.kilometer))


# In[22]:


GEO.rcParams["solver"] = "mg"
GEO.rcParams["initial.nonlinear.tolerance"] = 2e-2
GEO.rcParams["nonlinear.tolerance"] = 2e-2

# In[23]:


Model.init_model()

import underworld.function as fn

def post_hook():
    coords = fn.input()
    zz = coords[0] / (GEO.nd(Model.maxCoord[0]) - GEO.nd(Model.minCoord[0]))
    fact = fn.math.pow(fn.math.tanh(zz*20.0) + fn.math.tanh((1.0-zz)*20.0) - fn.math.tanh(20.0), 4)
    Model.plasticStrain.data[:] = Model.plasticStrain.data[:] * fact.evaluate(Model.swarm)

Model.postSolveHook = post_hook

Model.run_for(nstep=10)
Model.checkpoint(0)

