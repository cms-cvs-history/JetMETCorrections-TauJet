#ifndef TCTauAlgorithm_H
#define TCTauAlgorithm_H

/** \class TCTauAlgo
 *
 * Calculates TCTau based on detector response to charged particles
 * using the tracker 
 *
 * \authors    R.Kinnunen, S.Lehti, A.Nikitenko
 *
 * \version   1st Version July 2, 2009
 ************************************************************/

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TLorentzVector.h"
#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"

#include "FWCore/Framework/interface/Event.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Math/VectorUtil.h"

using namespace ROOT::Math;
using namespace std;
using namespace reco;
using namespace edm;

class TCTauAlgorithm {

    public:
        TCTauAlgorithm();
	TCTauAlgorithm(const edm::ParameterSet&);
        ~TCTauAlgorithm();

	TLorentzVector recalculateEnergy(const reco::CaloTau&);
        TLorentzVector recalculateEnergy(const reco::CaloJet&);
        TLorentzVector recalculateEnergy(const reco::Jet&);
        TLorentzVector recalculateEnergy(const reco::IsolatedTauTagInfo&);
	TLorentzVector recalculateEnergy(const reco::CaloJet&,const TrackRef&,const TrackRefVector&);

	void inputConfig(const edm::ParameterSet& iConfig);
	void eventSetup(const edm::Event&,const edm::EventSetup&);

	double efficiency();
	int    allTauCandidates();
	int    statistics();

    private:
        void init();

	XYZVector                 trackEcalHitPoint(const TransientTrack&,const CaloJet&);
	pair<XYZVector,XYZVector> getClusterEnergy(const CaloJet&,XYZVector&,double);
	XYZVector                 getCellMomentum(const CaloCellGeometry*,double&);


        int all,
            passed;

	int prongs;

	const TransientTrackBuilder* transientTrackBuilder;

	double signalCone;
	double ecalCone;
	double matchingCone;
	double tkptmin;
	double tkmaxipt;
	double tkmaxChi2;
	int    tkminPixelHitsn;
	int    tkminTrackerHitsn;

	double  etCaloOverTrackMin,
		etCaloOverTrackMax,
		etHcalOverTrackMin,
		etHcalOverTrackMax;

        InputTag trackInput;
        InputTag vertexInput;

        InputTag EcalRecHitsEB_input;
        InputTag EcalRecHitsEE_input;
        InputTag HBHERecHits_input;
        InputTag HORecHits_input;
        InputTag HFRecHits_input;

        const CaloSubdetectorGeometry* EB;
        const CaloSubdetectorGeometry* EE;
        const CaloSubdetectorGeometry* HB;
        const CaloSubdetectorGeometry* HE;
        const CaloSubdetectorGeometry* HO;
        const CaloSubdetectorGeometry* HF;

        Handle<EBRecHitCollection>   EBRecHits;
        Handle<EERecHitCollection>   EERecHits;

        Handle<HBHERecHitCollection> HBHERecHits;
        Handle<HORecHitCollection>   HORecHits;
        Handle<HFRecHitCollection>   HFRecHits;

	Handle<TrackCollection> tracks;
	Handle<VertexCollection> thePVs;
};
#endif

