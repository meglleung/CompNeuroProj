# Adapted by Aman Aberra from NEURON python tutorial:
# https://github.com/neuronsimulator/nrn/blob/master/docs/tutorials/ballandstick.py
from neuron import h
import numpy as np
h.load_file("stdrun.hoc")

def range_assignment(sec, var, start, stop):
    """linearly assign values between start and stop to each segment.var in section"""
    delta = stop - start
    for seg in sec:
        setattr(seg, var, start + seg.x * delta)
        
class Cell:
    def __init__(self, gid,ais_mode="soma",acd_connect_x=0.1,ais_length=60):
        self._gid = gid
        self.ais_mode = ais_mode
        self.aid_length = ais_length
        self._setup_morphology(ais_mode=ais_mode,acd_connect_x=acd_connect_x,ais_length=ais_length)
        self.all = self.soma.wholetree()
        self._setup_biophysics()        
        # h.define_shape()    
        self._setup_recordings()
        self.num_segments = np.sum([sec.nseg for sec in self.all])
    def __repr__(self):
        return "{}[{}]".format(self.name, self._gid)
    

class BallAndStick(Cell):
    name = "BallAndStick"

    def _setup_morphology(self,ais_mode="soma",acd_connect_x=0.1,apic_diam0=2.5,apic_diam1=0.5,ais_length=60):  # 4 to 2 for Hodepp
        self.soma = h.Section(name="soma", cell=self)
        self.apic = []
        # Main apical branch
        self.apic.append(h.Section(name="apic[0]", cell=self))
        self.apic[0].connect(self.soma,0.5)            
        # Apical tuft
        self.apic.append(h.Section(name="apic[1]", cell=self))
        self.apic.append(h.Section(name="apic[2]", cell=self))
        self.apic[1].connect(self.apic[0],1) # connect apical tuft branches
        self.apic[2].connect(self.apic[0],1) # connect apical tuft branches
        # Basal dendrites
        self.dend = []
        self.dend.append(h.Section(name="dend[0]",cell=self))
        self.dend.append(h.Section(name="dend[1]",cell=self))
        self.dend[0].connect(self.soma,1)
        self.dend[1].connect(self.soma,1)

        # soma geometry
        # self.soma.pt3dadd(-6,0,0,8) # L = 12, diam = 8 - Hodapp 2022
        # self.soma.pt3dadd(6,0,0,8)
        self.soma.pt3dadd(-10,0,0,20) # L = 20, diam = 20 - Gulledge and Bravo 2016
        self.soma.pt3dadd(10,0,0,20)
        self.soma.nseg = 5
        # apical dendrite geometry                 
        self.apic[0].pt3dadd(0,0,0,apic_diam0)
        self.apic[0].pt3dadd(0,400,0,apic_diam1) 
        self.apic[0].nseg = 1 + 2*int(self.apic[0].L/40) # 1 segments per 2 µm 
        # for seg in self.apic[0]:
        #     # taper from 4 to 2 µm
        #     seg.diam = np.interp(seg.x, [0, 1], [apic_diam0, apic_diam1]) 
        # apical tuft geometry
        self.apic[1].L = 150
        self.apic[1].diam = 2
        self.apic[1].nseg = 5
        self.apic[2].L = 150
        self.apic[2].diam = 2        
        self.apic[2].nseg = 5
        # basal dendrite geometry
        self.dend[0].pt3dadd(-6,0,0,1.8) # L = 200, diam = 1.8
        self.dend[0].pt3dadd(-206,0,0,1.8)
        self.dend[0].nseg = 1 + 2*int(self.dend[0].L/2) # 1 segments per 2 µm 
        self.dend[1].pt3dadd(6,0,0,1.8)
        self.dend[1].pt3dadd(206,0,0,1.8)     
        self.dend[1].nseg = 1 + 2*int(self.dend[0].L/40) # 1 segments per 40 µm 
        # Axon geometry -- 2 sections so would have to be /2
        self.ais_prox = h.Section(name="ais_prox",cell=self)
        self.ais_prox.L = ais_length/2 # 17 - Hodapp
        self.ais_prox.diam =  1.5 # 1.22 Hodapp
        self.ais_prox.nseg = 21       
        self.ais_dist = h.Section(name="ais_dist",cell=self)        
        self.ais_dist.connect(self.ais_prox)
        self.ais_dist.L = ais_length/2 # 17 - Hodapp
        self.ais_dist.diam =  1.5 # 1.22 Hodapp
        self.ais_dist.nseg = 21
        self.axon = h.Section(name="axon",cell=self)        
        self.axon.connect(self.ais_dist)
        self.axon.diam = 1
        self.axon.L = 500
        self.axon.nseg = 1 + 2*int(self.axon.L/40) # 1 segments per 40 µm 
        if ais_mode == 'soma': # connect to soma
            self.ais_prox.connect(self.soma)
        elif ais_mode == 'dend': # connect to basal dendritic branch
            self.ais_prox.connect(self.dend[0],acd_connect_x)
        else:
            raise AssertionError('ais_mode should be \'soma\' or \'dend\'')
        

    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100  #  Ohm * cm 100 - Gulledge, 200 - Hodapp
            sec.cm = 1  # micro Farads / cm^2
            sec.insert('extracellular') # allows recording membrane current
            # Make both soma and dendrite passive with same parameters
            sec.insert('pas')
            sec.g_pas = 1/15e3 # S/cm2 - Rm = 25e3 Ohm*cm2 - Hodapp 
            sec.e_pas = -70 # mV     
            # active conductances
            sec.insert('na')                    
            sec.insert('kv')
        # Set somatic conductances
        self.soma.gbar_na = 400 # 100 pS/um2
        self.soma.gbar_kv = 100

        # soma_mechs = ['inaT','ikdT','imZ','hNa','kap','cal','cat','somacar','kca','mykca','cad']
        # [self.soma.insert(m) for m in soma_mechs] # insert mechanisms
        # self.soma.gnabar_inaT = 0.2 # mS/cm2
        # self.soma.vtraub_inaT = -58 # 

        # # Proximal AIS
        # ais_prox_mechs = ['inaT','ikdT','imZ']
        # [self.ais_prox.insert(m) for m in ais_prox_mechs]
        self.ais_prox.gbar_na = 100
        self.ais_prox.gbar_kv = 100
        # # Distal AIS
        # ais_dist_mechs = ['inaT','ikdT','imZ'] 
        # [self.ais_dist.insert(m) for m in ais_dist_mechs]
        self.ais_dist.gbar_na = 8000
        self.ais_dist.gbar_kv = 2000
        # # Main axon
        # axon_mechs = ['inaT','ikdT']        
        # [self.axon.insert(m) for m in axon_mechs]
        self.axon.gbar_na=300
        self.axon.gbar_kv=60
        # main apical dendritic branch
        # apic_mechs = ['hha_old','imZ','hNa','kap','calH','cat','car','kca','mykca','cad']
        # [self.apic[0].insert(m) for m in apic_mechs]
        range_assignment(self.apic[0],'gbar_na',100,20)
        range_assignment(self.apic[0],'gbar_kv',100,20)        
        self.apic[1].gbar_na = 20
        self.apic[2].gbar_na = 20
        self.apic[1].gbar_kv = 20
        self.apic[2].gbar_kv = 20
        # apical tuft 
        # apic_tuft_mechs = ['hNa','kad','hha_old']
        # [sec.insert(m) for m in apic_tuft_mechs for sec in self.apic[1:]]

        # Basal dendritic branches
        range_assignment(self.dend[0],'gbar_na',100,20)
        range_assignment(self.dend[0],'gbar_kv',100,20)
        range_assignment(self.dend[1],'gbar_na',100,20)
        range_assignment(self.dend[1],'gbar_kv',100,20)        
        # dend_mechs = ['cad','mykca','calH','car','cat','hNa','hha_old','kap','kca']
        # [sec.insert(m) for m in dend_mechs for sec in self.dend]
        
        for sec in self.all:
            sec.ena = 50
            sec.ek = -80

    def _setup_recordings(self):
        self._spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma) # monitor V in soma for spiking
        self.spike_times = h.Vector() 
        self._spike_detector.record(self.spike_times) # Record spike times to spike_times Vector        
        self.soma_v = h.Vector().record(self.soma(0.5)._ref_v) # Record somatic voltage        
        self.ais_v = h.Vector().record(self.ais_dist(0.5)._ref_v) # Record from distal AIS
        self.dend0_vs = [h.Vector().record(seg._ref_v) for seg in self.dend[0]] # Record voltages from basal dendrite carrying axon
        self.dend1_vs = [h.Vector().record(seg._ref_v) for seg in self.dend[1]] # Record voltages from other basal dendrite
