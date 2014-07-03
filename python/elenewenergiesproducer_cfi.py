## 
#This module is a producer of ValueMaps that associate to each electron
#in the GsfElectron collection, one float that is a new type of energy
#associated to the electron. 
#The core of the module is the plugin EleNewEnergiesProducer.cc whose
#cfi.py file is elenewenergiesproducer_cfi.py
#
#The plugin just create the new collection, loop over the gsfElectron
#collection and calculate the new energy according to a separate class
#that should called by the plugin. The implementation of the
#calculation of the new energy should not be done in the plugin but in
#a separate class in order to keep it more general and flexible.
#The separate class can be added in the src directory or already
#present in some other CMSSW package (add the required lines in the
#BuildFile.xml).
#
#
#You should copy the files indicated in data/copy.url in the data/ directory
##
import FWCore.ParameterSet.Config as cms

# ALCARECO collections
eleNewEnergiesProducer = cms.EDProducer('EleNewEnergiesProducer',
                                        electronCollection = cms.InputTag('patPhotons'),
                                        recHitCollectionEB = cms.InputTag("reducedEcalRecHitsEB"),
                                        recHitCollectionEE = cms.InputTag("reducedEcalRecHitsEE"),
                                        rhoFastJet = cms.InputTag('kt6PFJets',"rho"),
                                        vertexCollection = cms.InputTag('offlinePrimaryVerticesWithBS'),
                                        BeamSpotCollection = cms.InputTag('offlineBeamSpot'),
                                        conversionCollection = cms.InputTag('patConversionsPFlowAK5chs'),
                                        isMC = cms.bool(False),
                                        regrPhoFile = cms.string('./data/gbrv3ph_52x.root'),
                                        regrEleFile = cms.string('./data/gbrv3ele_52x.root'),
                                        regrEleFile_fra = cms.string('nocorrection.root'),
                                        correctionFileName = cms.string(''),
                                        correctionType = cms.string(''),
                                        regrPhoJoshV5_SemiParamFile = cms.string('src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v5_forest_ph.root'),
                                        regrEleJoshV5_SemiParamFile = cms.string('src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v5_forest_ele.root'),

                                        )
