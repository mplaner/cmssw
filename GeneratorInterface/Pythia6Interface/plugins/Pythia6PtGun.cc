
#include <iostream>

#include "Pythia6PtGun.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "GeneratorInterface/Pythia6Interface/interface/PYR.h"

#include "FWCore/Framework/interface/MakerMacros.h"

using namespace edm;
using namespace gen;

Pythia6PtGun::Pythia6PtGun( const ParameterSet& pset ) :
   Pythia6Gun(pset)
{
   
   // ParameterSet defpset ;
   ParameterSet pgun_params = 
      pset.getParameter<ParameterSet>("PGunParameters"); //, defpset ) ;
   fMinEta     = pgun_params.getParameter<double>("MinEta"); // ,-2.2);
   fMaxEta     = pgun_params.getParameter<double>("MaxEta"); // , 2.2);
   fMinPt       = pgun_params.getParameter<double>("MinPt"); // ,  20.);
   fMaxPt       = pgun_params.getParameter<double>("MaxPt"); // , 420.);
   fAddAntiParticle = pgun_params.getParameter<bool>("AddAntiParticle"); //, false) ;  

}

Pythia6PtGun::~Pythia6PtGun()
{
}

void Pythia6PtGun::generateEvent()
{
   
   Pythia6Service::InstanceWrapper guard(fPy6Service);	// grab Py6 instance

   // now actualy, start cooking up the event gun 
   //

   // 1st, primary vertex
   //
   HepMC::GenVertex* Vtx = new HepMC::GenVertex( HepMC::FourVector(0.,0.,0.));

   // here re-create fEvt (memory)
   //
   fEvt = new HepMC::GenEvent() ;
     
   int ip=1;
   for ( size_t i=0; i<fPartIDs.size(); i++ )
   {
	 int particleID = fPartIDs[i]; // this is PDG - need to convert to Py6 !!!
	 int py6PID = HepPID::translatePDTtoPythia( particleID );
         int dum = 0;
	 double pt=0, mom=0, ee=0, the=0, eta=0;
	 double mass = pymass_(particleID);
	 double phi = (fMaxPhi-fMinPhi)*pyr_(&dum)+fMinPhi;

	 eta  = (fMaxEta-fMinEta)*pyr_(&dum)+fMinEta;                                                      

	 the  = 2.*atan(exp(-eta));                                                                          
	 
	 pt = (fMaxPt-fMinPt)*pyr_(&dum)+fMinPt;
	 
	 mom = pt/sin(the);
	 ee = sqrt(mom*mom+mass*mass);

	 py1ent_(ip, py6PID, ee, the, phi);
	 
         double px     = pt*cos(phi) ;
         double py     = pt*sin(phi) ;
         double pz     = mom*cos(the) ;
         
	 HepMC::FourVector p(px,py,pz,ee) ;
         HepMC::GenParticle* Part = 
             new HepMC::GenParticle(p,particleID,1);
         Part->suggest_barcode( ip ) ;
         Vtx->add_particle_out(Part);
	 
	 if(fAddAntiParticle)
	 {
	    ip = ip + 1;
// Check if particle is its own anti-particle.
            int pythiaCode = pycomp_(py6PID); // this is py6 internal validity check, it takes Pythia6 pid
	                                          // so actually I'll need to convert
            int has_antipart = pydat2.kchg[3-1][pythiaCode-1];
            int particleID2 = has_antipart ? -1 * particleID : particleID;
	    int py6PID2 = has_antipart ? -1 * py6PID : py6PID;	 // this py6 id, for py1ent    
	    the = 2.*atan(exp(eta));
	    phi  = phi + M_PI;
	    if (phi > 2.* M_PI) {phi = phi - 2.* M_PI;}         
	    py1ent_(ip, py6PID2, ee, the, phi);
            HepMC::FourVector ap(-px,-py,-pz,ee) ;
	    HepMC::GenParticle* APart =
	       new HepMC::GenParticle(ap,particleID2,1);
	    APart->suggest_barcode( ip ) ;
	    Vtx->add_particle_out(APart) ;	    
	 }
	 ip++;
   }
   
   fEvt->add_vertex(Vtx);
     
   // run pythia
   pyexec_();
   
   return;
}

DEFINE_FWK_MODULE(Pythia6PtGun);
