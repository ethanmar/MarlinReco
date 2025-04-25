#include "TauFinder.h"
#include <iostream>
#include <iomanip>

using namespace std;


#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCRelation.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include "UTIL/LCRelationNavigator.h"
#include <EVENT/Vertex.h>

#include <marlin/Global.h>

#include "HelixClass.h"
#include "SimpleLine.h"
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#define coutEv -1

#define coutUpToEv 0


using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

TauFinder aTauFinder ;


bool MyEnergySort( ReconstructedParticle *p1, ReconstructedParticle *p2)
{
  return fabs(p1->getEnergy()) > fabs(p2->getEnergy());
}

// Function to calculate pseudorapidity (eta)
double computeEta(const double p[3]) {
  double pz = p[2];
  double pt = sqrt(p[0] * p[0] + p[1] * p[1]);
  return 0.5 *
         log((sqrt(pt * pt + pz * pz) + pz) / (sqrt(pt * pt + pz * pz) - pz));
}

// Function to calculate azimuthal angle (phi)
double computePhi(const double p[3]) { return atan2(p[1], p[0]); }

// Function to compute deltaR between two particles
double computeDeltaR(const double pvec[3], const double pvec_tau[3]) {
  double eta1 = computeEta(pvec);
  double eta2 = computeEta(pvec_tau);

  double phi1 = computePhi(pvec);
  double phi2 = computePhi(pvec_tau);

  double deltaEta = eta1 - eta2;
  double deltaPhi = phi1 - phi2;

  // Ensure deltaPhi is within the range [-π, π]
  if (deltaPhi > M_PI)
    deltaPhi -= 2 * M_PI;
  if (deltaPhi < -M_PI)
    deltaPhi += 2 * M_PI;

  return sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
}

TauFinder::TauFinder() : Processor("TauFinder") 
{
  // modify processor description
  _description = "TauFinder writes tau candidates as ReconstructedParticles into collection. It runs on a collection of ReconstructedParticels, if you want  to run on MCParticles you have to convert them before hand (use e.g. PrepareRECParticles processor)" ;
  
  // register steering parameters: name, description, class-variable, default value
  
 
  //the link between the reconstructed tau and the input particles used for it
 

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "PFOCollection",
                            "Collection of PFOs",
                            _incol ,
                            std::string("PandoraPFANewPFOs"));

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "TauRecCollection",
                            "Collection of Tau Candidates",
                            _outcol ,
                            std::string("TauRec_PFO"));
  
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "TauRecRestCollection",
                            "Collection of Particles in Rest Group not in Tau Candidates",
                            _outcolRest ,
                            std::string("TauRecRest_PFO"));
  
  registerOutputCollection( LCIO::LCRELATION,
			    "TauRecLinkCollectionName" , 
			    "Name of the Tau link to ReconstructedParticle collection"  ,
			    _colNameTauRecLink ,
			    std::string("TauRecLink_PFO") ) ;
  
  registerProcessorParameter( "FileName_Signal",
			      "Name of the Signal output file "  ,
			      _OutputFile_Signal ,
			      std::string("Signal.root") ) ;

 registerProcessorParameter( "pt_cut" ,
                              "Cut on pt to suppress background"  ,
                              _ptcut ,
                              (float)0.2) ;
  registerProcessorParameter( "cosT_cut" ,
                              "Cut on cosT to suppress background"  ,
                              _cosTcut ,
                              (float)0.99) ;
  
  registerProcessorParameter( "searchConeAngle" ,
                              "Opening angle of the search cone for tau jet in rad"  ,
                              _coneAngle ,
                              (float)0.05) ;

  registerProcessorParameter( "isolationConeAngle" ,
                              "Outer isolation cone around search cone of tau jet in rad (relativ to cone angle)"  ,
                              _isoAngle ,
                              (float)0.02) ;
  
  registerProcessorParameter( "isolationEnergy" ,
                              "Energy allowed within isolation cone region"  ,
                              _isoE ,
                              (float)5.0) ;
  
  registerProcessorParameter( "ptseed" ,
                              "Minimum tranverse momentum of tau seed"  ,
                              _ptseed ,
                              (float)5.0) ;
  
  registerProcessorParameter( "invariant_mass" ,
                              "Upper limit on invariant mass of tau candidate"  ,
                              _minv ,
                              (float)2.0) ;
}


void TauFinder::init() 
{ 
  streamlog_out(DEBUG) << " init called  " << std::endl ;
  
  // usually a good idea to
  printParameters() ;

  rootfile = new TFile((_OutputFile_Signal).c_str (),"RECREATE");

  anatree = new TTree("anatree", "TTree with information on taufinder process");

  if (anatree == 0) {
    throw lcio::Exception(" Invalid tree pointer !!! ");
  }
  anatree->Branch(("ntau"), &_ntau, ("ntau/I"));

  anatree->Branch(("ngood"), &_ngood, ("ngood/I"));
  anatree->Branch(("nfailseed"), &_nfail_seed, ("nfailseed/I"));
  anatree->Branch(("nfailQtrack"), &_nfail_Qtrack, ("nfailQtrack/I"));
  anatree->Branch(("nrej"), &_nrej, ("nrej/I"));

  anatree->Branch(TString("t_isoE"), _tau_isoE, "t_isoE[ntau]/F");
  anatree->Branch(TString("t_minv"), _tau_minv, TString("t_minv[ntau]/F"));
  anatree->Branch(TString("t_pt"), _tau_pt, TString("t_pt[ntau]/F"));
  anatree->Branch(TString("t_p"), _tau_p, TString("t_p[ntau]/F"));
  anatree->Branch(TString("t_ene"), _tau_ene, TString("t_ene[ntau]/F"));
  anatree->Branch(TString("t_phi"), _tau_phi, TString("t_phi[ntau]/F"));
  anatree->Branch(TString("t_eta"), _tau_eta, TString("t_eta[ntau]/F"));
  anatree->Branch(TString("t_nQ"), _tau_nQ, TString("t_nQ[ntau]/F"));
  anatree->Branch(TString("t_nN"), _tau_nN, TString("t_nN[ntau]/F"));
  anatree->Branch(TString("t_nQN"), _tau_nQN, TString("t_nQN[ntau]/F"));
  anatree->Branch(TString("t_good"), _tau_good, TString("t_good[ntau]/I"));

  anatree->Branch(("nrej_isoE"), &_nrej_isoE, ("nrej_isoE/I"));
  anatree->Branch(("nrej_minv"), &_nrej_minv, ("nrej_minv/I"));
  anatree->Branch(("nrej_nQ"), &_nrej_nQ, ("nrej_nQ/I"));
  anatree->Branch(("nrej_nQN"), &_nrej_nQN, ("nrej_nQN/I"));

  anatree->Branch(("event_num"), _event_num, ("event_num[ntau]/I"));

  _nRun = 0 ;
  _nEvt = 0 ;
  _mergeTries=0;
}

void TauFinder::processRunHeader( LCRunHeader* )
{ 
  _nRun++ ;
} 

void TauFinder::processEvent( LCEvent * evt ) 
{ 

  // this gets called for every event 
  // usually the working horse ...

  _nrej = 0;
  _nrej_isoE = 0;
  _nrej_minv = 0;
  _nrej_nQ = 0;
  _nrej_nQN = 0;
  _ntau = 0;
  _ngood = 0;
  _nfail_seed = 0;
  _nfail_Qtrack = 0;
  
  LCCollection *colRECO;

  try {
    colRECO = evt->getCollection( _incol ) ;
} catch (Exception& e) {
    colRECO = 0;
  }
 
  //LCRelation: to keep information from which particles the tau was made
  LCCollectionVec *relationcol = new LCCollectionVec(LCIO::LCRELATION);
  relationcol->parameters().setValue(std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  relationcol->parameters().setValue(std::string("ToType"),LCIO::RECONSTRUCTEDPARTICLE);
 
  _nEvt = evt->getEventNumber();  
  
  LCCollectionVec * reccol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec * restcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  restcol->setSubset(1);
  //sort all input particles into charged and neutral
  std::vector<ReconstructedParticle*> Avector;//all particles
  std::vector<ReconstructedParticle*> Qvector;//charged particles
  std::vector<ReconstructedParticle*> Nvector;//neutral particles


  if( colRECO != 0)
    {
      int nRCP = colRECO->getNumberOfElements();
  
      for (int i = 0; i < nRCP; i++) 
	{
	  ReconstructedParticle *particle = static_cast<ReconstructedParticle*>( colRECO->getElementAt( i ) );
	  double pt=sqrt(particle->getMomentum()[0]*particle->getMomentum()[0]
			 +particle->getMomentum()[1]*particle->getMomentum()[1]);
	  double Cos_T  = fabs(particle->getMomentum()[2]) / sqrt(pow(particle->getMomentum()[0],2)+pow(particle->getMomentum()[1],2) + pow(particle->getMomentum()[2],2));
	  if(pt<_ptcut || Cos_T>_cosTcut)   
	    continue;
	  Avector.push_back(particle);
	  if(particle->getCharge()==1 || particle->getCharge()==-1)
	    Qvector.push_back(particle);
	  else
	    Nvector.push_back(particle);
	}
    }//colRECO

  
  //sort mc vec according to energy
  std::sort(Qvector.begin(), Qvector.end(), MyEnergySort);
  std::sort(Nvector.begin(), Nvector.end(), MyEnergySort);
  
  //vector to hold tau candidates
  std::vector<std::vector<ReconstructedParticle*> > tauvec;
  bool finding_done=false;
  while(Qvector.size() && !finding_done)
    finding_done= FindTau(Qvector,Nvector,tauvec);

  if (tauvec.size() == 0) {

    if (Qvector.size() == 0) {
      if (_nEvt < coutUpToEv || _nEvt == coutEv)
        streamlog_out(DEBUG)
            << "No charged particle found in Evt " << _nEvt << endl;
      ++_nfail_Qtrack;
    }

    else
      ++_nfail_seed;
  }
  
  //combine associated particles to tau
  std::vector<std::vector<ReconstructedParticle*> >::iterator iterT=tauvec.begin();
  std::vector<ReconstructedParticleImpl* > tauRecvec;
  //remember number of charged tracks in each tau for possible merging later
  std::vector<int> QTvec;
  std::vector<int> NTvec;
  for(unsigned int p=0;p<tauvec.size();p++)
    {
      ReconstructedParticleImpl *taurec=new ReconstructedParticleImpl();
      double E=0;
      double charge=0;
      double mom[3]={0,0,0};
      int chargedtracks=0;
      int neutraltracks=0;
      std::vector<ReconstructedParticle*> tau=tauvec[p];
     
      for(unsigned int tp=0;tp<tau.size();tp++)
	{
	  //add up energy and momentum
	  E+=tau[tp]->getEnergy();
	  charge+=tau[tp]->getCharge();
	  mom[0]+=tau[tp]->getMomentum()[0];
	  mom[1]+=tau[tp]->getMomentum()[1];
	  mom[2]+=tau[tp]->getMomentum()[2];
	  if(tau[tp]->getCharge())
	    chargedtracks++;
	  else
	    neutraltracks++;
	  taurec->addParticle(tau[tp]);
	}
      
      double pt_tau=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
      double psquare=pt_tau*pt_tau+mom[2]*mom[2];
      double mass_inv=0;
      
      if(E*E<psquare)
	mass_inv=E-sqrt(psquare);
      else
	mass_inv=sqrt(E*E-psquare);
     
      //check for invariant mass
      if(mass_inv>_minv || mass_inv<-0.001 || chargedtracks>4 || chargedtracks==0)
	{
	  if(mass_inv>_minv){
	    _nrej_minv++;
	  }
	  if(mass_inv<-0.001){
	    _nrej_minv++;
	  }
	  if(chargedtracks>4 || chargedtracks==0){
	    _nrej_nQ++;
	  }

	  //put those particles into rest group
	  for(unsigned int tp=0;tp<tau.size();tp++)
	    restcol->addElement(tau[tp]);
	  
	  iterT=tauvec.erase(iterT);
	  p--;
	  
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    {
	      double phi=180./TMath::Pi()*atan(mom[1]/mom[0]);
	      double theta=180./TMath::Pi()*atan(pt_tau/fabs(mom[2])); 
	      streamlog_out(DEBUG) <<"Tau candidate failed: minv="<<mass_inv<<"   pt="<<pt_tau<<" "<<E<<" Q trks:"<<chargedtracks<<" N trks:"<<neutraltracks<<" "<<phi<<" "<<theta<<endl;
	    }
	  delete taurec;
	  continue;
	}
      else
	++iterT;

      for(unsigned int tp=0;tp<tau.size();tp++) {
	LCRelationImpl *rel = new LCRelationImpl(taurec,tau[tp]);
	relationcol->addElement( rel );
      }
      
      int pdg=15;
      if(charge<0)
	pdg=-15;
	
      taurec->setEnergy(E);
      taurec->setCharge(charge);
      taurec->setMomentum(mom);
      taurec->setType(pdg);
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	{
	  double phi=180./TMath::Pi()*atan(mom[1]/mom[0]);
	  double theta=180./TMath::Pi()*atan(pt_tau/fabs(mom[2])); 
	  streamlog_out(DEBUG)<<"Tau candidate "<<p<<": "<<E<<" Q trks:"<<chargedtracks<<" N trks:"<<neutraltracks<<" "<<phi<<" "<<theta<<endl;
	}
      QTvec.push_back(chargedtracks);
      NTvec.push_back(neutraltracks);
      tauRecvec.push_back(taurec);
    }
  //merge taus that are very close together, because they are likely to be from 1 tau that got split in algorithm
  LCRelationNavigator *relationNavigator = new LCRelationNavigator( relationcol );

  if(tauRecvec.size()>1)
    {
      std::vector<ReconstructedParticleImpl*>::iterator iterC=tauRecvec.begin();
      std::vector<ReconstructedParticleImpl*>::iterator iterF=tauRecvec.begin();
      int erasecount=0;
      for ( unsigned int t=0; t<tauRecvec.size() ; t++ )
	{
	  ReconstructedParticleImpl *tau=static_cast<ReconstructedParticleImpl*>(tauRecvec[t]);
	  const double *mom=tau->getMomentum();
	  double pt_tau=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
	  double phi=180./TMath::Pi()*atan(mom[1]/mom[0]);
	  double theta=180./TMath::Pi()*atan(pt_tau/fabs(mom[2])); 
	  double angle=10000000;
	 
	  for ( unsigned int t2=t+1; t2<tauRecvec.size() ; t2++ )
	    {
	      iterC=tauRecvec.begin()+t2;
	      ReconstructedParticleImpl *taun=static_cast<ReconstructedParticleImpl*>(tauRecvec[t2]);
	      
	      const double *momn=taun->getMomentum();
	      angle=acos((mom[0]*momn[0]+mom[1]*momn[1]+mom[2]*momn[2])/
			 (sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])*
			  sqrt(momn[0]*momn[0]+momn[1]*momn[1]+momn[2]*momn[2])));
	      //merging happens
	      if(angle<_coneAngle)
		{
		  _mergeTries++;
		  double E=tau->getEnergy();
		  double En=E+taun->getEnergy();
		  tau->setEnergy(En);
		  double newp[3]={mom[0]+momn[0],mom[1]+momn[1],mom[2]+momn[2]};
		  tau->setMomentum(newp);		  
		  tau->setCharge(tau->getCharge()+taun->getCharge());

		  if(_nEvt<coutUpToEv || _nEvt==coutEv)
		    {
		      streamlog_out(DEBUG)<<" Tau Merging: "<<endl;
		      streamlog_out(DEBUG)<<t<<" "<<E<<" "<<phi<<" "<<theta<<endl;
		      streamlog_out(DEBUG)<<t2<<" "<<taun->getEnergy()<<" "<<angle<<" -> "<<En<<" "<<QTvec[t]+QTvec[t2]<<endl;
		    }
		  
		  //check for invariant mass and number of tracks
		  double ptn=sqrt(newp[0]*newp[0]+newp[1]*newp[1]);
		  double psquaren=ptn*ptn+newp[2]*newp[2];
		  double mass_inv=0;
		  if(En*En<psquaren)
		    mass_inv=En*En-psquaren;
		  else
		    mass_inv= sqrt(En*En-psquaren);
		  //failed to merge
		  if(mass_inv>_minv || mass_inv<-0.001 ||  QTvec[t+erasecount]+QTvec[t2+erasecount]>4)
		    {
		      if(mass_inv>_minv){
			_nrej_minv++;
		      }
		      if(mass_inv<-0.001){
			_nrej_minv++;
		      }
		      if(QTvec[t+erasecount]+QTvec[t2+erasecount]>4){
       			_nrej_nQ++;
		      }
		      
		      //put those particles into rest group
		      for(unsigned int tp=0;tp<(*iterC)->getParticles().size();tp++)
			restcol->addElement((*iterC)->getParticles()[tp]);
		      for(unsigned int tp=0;tp<(*iterF)->getParticles().size();tp++)
			restcol->addElement((*iterF)->getParticles()[tp]);
		      
		      delete *iterC;
		      tauRecvec.erase(iterC);
		      erasecount++;
		      t2--;
		      delete *iterF;
		      iterF=tauRecvec.erase(iterF);
		      erasecount++;
		      if(iterF!=tauRecvec.end())
			{
			  tau=*iterF;
			  mom=tau->getMomentum();
			  pt_tau=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
			  phi=180./TMath::Pi()*atan(mom[1]/mom[0]);
			  theta=180./TMath::Pi()*atan(pt_tau/fabs(mom[2])); 
			  angle=10000000;
			}
		    }
		  else //merge
		    {
		      //set the relations and add particles from one tau to the other
		      std::vector< ReconstructedParticle * > mergetaus=taun->getParticles();
		      for(unsigned int p=0;p<mergetaus.size();p++)
			tau->addParticle(mergetaus[p]);
		      EVENT::LCObjectVec relobjFROM = relationNavigator->getRelatedToObjects(taun);
		      for(unsigned int o=0;o<relobjFROM.size();o++)
			{
			  ReconstructedParticle *rec=static_cast<ReconstructedParticle*>(relobjFROM[o]);
			  LCRelationImpl *rel = new LCRelationImpl(tau,rec);
			  relationcol->addElement( rel );
			}
		      delete *iterC;
		      tauRecvec.erase(iterC);
		      erasecount++;
		      t2--;
		    }
		}
	    }
	  iterF++;
	}
    }
  delete relationNavigator;
  //test for isolation and too many tracks
  std::vector<ReconstructedParticleImpl*>::iterator iter=tauRecvec.begin();
  _ntau = tauRecvec.size();
  int erasecount=0;
  for ( unsigned int t=0; t<tauRecvec.size() ; t++ )
    {
      ReconstructedParticleImpl *tau=static_cast<ReconstructedParticleImpl*>(tauRecvec[t]);
      double E_iso=0;
      int nparticles=0;
      const double *pvec_tau=tau->getMomentum();
      const double p_tau =
        sqrt(pvec_tau[0] * pvec_tau[0] + pvec_tau[1] * pvec_tau[1] +
             pvec_tau[2] * pvec_tau[2]);

      _tau_nQ[t + erasecount] = QTvec[t + erasecount];
      _tau_nN[t + erasecount] = NTvec[t + erasecount];
      _tau_nQN[t + erasecount] = QTvec[t + erasecount] + NTvec[t + erasecount];
      _tau_ene[t + erasecount] = tau->getEnergy();
      _tau_p[t + erasecount] = p_tau;
      _tau_pt[t + erasecount] =
        sqrt(pvec_tau[0] * pvec_tau[0] + pvec_tau[1] * pvec_tau[1]);
      _tau_phi[t + erasecount] = TMath::ATan2(pvec_tau[1], pvec_tau[0]);
      _tau_eta[t + erasecount] =
        0.5 * TMath::Log((p_tau + pvec_tau[2]) / (p_tau - pvec_tau[2]));
      ;
      _tau_minv[t + erasecount] =
        sqrt(tau->getEnergy() * tau->getEnergy() - p_tau * p_tau);
      
      //too many particles in tau 
      if(QTvec[t+erasecount]+NTvec[t+erasecount]>10 || QTvec[t+erasecount]>4)
	{
       	  _nrej_nQ++;
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    streamlog_out(DEBUG) <<"Tau "<<tau->getEnergy()
                                 <<": too many particles: "<<QTvec[t+erasecount]<<" "<<NTvec[t+erasecount]<<endl;
	  _tau_good[t + erasecount] = 0;
	  //put those particles into rest group
	  for(unsigned int tp=0;tp<(*iter)->getParticles().size();tp++)
	    restcol->addElement((*iter)->getParticles()[tp]);
	  delete *iter;
	  iter=tauRecvec.erase(iter);
	  erasecount++;
	  t--;
	  continue;
	}
      //isolation
      for ( unsigned int s=0; s<Avector.size() ; s++ )
	{
	  ReconstructedParticle *track=static_cast<ReconstructedParticle*>(Avector[s]);
	  const double *pvec=track->getMomentum();
	  double angle=acos((pvec[0]*pvec_tau[0]+pvec[1]*pvec_tau[1]+pvec[2]*pvec_tau[2])/
			    (sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2])*
			     sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]+pvec_tau[2]*pvec_tau[2])));
	  if(angle>_coneAngle && angle<_isoAngle+_coneAngle)
	    {
	      nparticles++;
	      E_iso+=track->getEnergy();
	    }
	}

      _tau_isoE[t + erasecount] = E_iso;

      _event_num[t + erasecount] = _nEvt;
      
      if(E_iso>_isoE)
	{
     	  _nrej_isoE++;
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    streamlog_out(DEBUG) <<"Tau "<<tau->getEnergy()
                                 <<": Isolation Energy: "<<E_iso<<" in "<<nparticles<<" particles"<<endl;
	  _tau_good[t+erasecount] = 0;
	  //put those particles into rest group
	  for(unsigned int tp=0;tp<(*iter)->getParticles().size();tp++)
	    restcol->addElement((*iter)->getParticles()[tp]);
	  delete *iter;
	  iter=tauRecvec.erase(iter);
	  erasecount++;
	  t--;
	}
      else
	{
	  reccol->addElement(tau);  
	  streamlog_out(DEBUG) << "Tau "<<tau->getEnergy()<<" "<<QTvec[t+erasecount]<<" "<<NTvec[t+erasecount]<<endl;	  
	  _tau_good[t+erasecount] = 1;
	  _ngood++;
	  iter++;
	}
    }

    _nrej = _nrej_isoE + _nrej_minv + _nrej_nQ + _nrej_nQN;
  
  //put remaining particles into rest group
  for(unsigned int tp=0;tp<Qvector.size();tp++)
    restcol->addElement(Qvector[tp]);
  for(unsigned int tp=0;tp<Nvector.size();tp++)
    restcol->addElement(Nvector[tp]);

  evt->addCollection(reccol,_outcol);
  evt->addCollection(restcol,_outcolRest);
  evt->addCollection(relationcol,_colNameTauRecLink);

  anatree->Fill();
  
  _nEvt ++ ;
  
}

bool TauFinder::FindTau(std::vector<ReconstructedParticle*> &Qvec,std::vector<ReconstructedParticle*> &Nvec,
			std::vector<std::vector<ReconstructedParticle*> > &tauvec)
{
  std::vector<ReconstructedParticle*>  tau;
  if(Qvec.size()==0)
    {
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	streamlog_out(DEBUG) << "No charged particle in event!"<<endl;
      return true;
    }
  double OpAngleMax=0;
  //find a good tauseed, check impact parameter 
  ReconstructedParticle *tauseed=NULL;
  std::vector<ReconstructedParticle*>::iterator iterS=Qvec.begin();
  for ( unsigned int s=0; s<Qvec.size() ; s++ )
    {
      tauseed=static_cast<ReconstructedParticle*>(Qvec[s]);
      float mom[3];
      for (int icomp=0; icomp<3; ++icomp) 
	mom[icomp]=(float)tauseed->getMomentum()[icomp];

      double pt=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);

      if(pt>_ptseed)
	break;
      else
	{
	  iterS++;
	  tauseed=NULL;
	}
    }
  if(!tauseed)
    {
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	streamlog_out(DEBUG) << "no further tau seed!"<<endl;
      return true;
    }
  
  double  Etau=tauseed->getEnergy();
 
  tau.push_back(tauseed);

  // just for printing out info
  {
    const double *pvec=tauseed->getMomentum();
    double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
    double p=sqrt(pt*pt+pvec[2]*pvec[2]);
    double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
    double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2]));
 
    if(_nEvt<coutUpToEv || _nEvt==coutEv)
      streamlog_out(DEBUG) << "seeding: "<<tauseed->getType()<<"\t"<<tauseed->getEnergy()<<"\t"<<p<<"\t"<<theta<<"\t"<<phi<<endl;
  }

  Qvec.erase(iterS);
  double pvec_tau[3]={0,0,0};
  pvec_tau[0]=tauseed->getMomentum()[0];
  pvec_tau[1]=tauseed->getMomentum()[1];
  pvec_tau[2]=tauseed->getMomentum()[2];

  //assign charged particles to tau candidate
  std::vector<ReconstructedParticle*>::iterator iterQ=Qvec.begin();
  for (unsigned int s=0; s<Qvec.size() ; s++ )
    {
      ReconstructedParticle *track=static_cast<ReconstructedParticle*>(Qvec[s]);

      const double *pvec=track->getMomentum();
      double angle=acos((pvec[0]*pvec_tau[0]+pvec[1]*pvec_tau[1]+pvec[2]*pvec_tau[2])/
			(sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2])*
			 sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]+pvec_tau[2]*pvec_tau[2])));
      double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
      double p=sqrt(pt*pt+pvec[2]*pvec[2]);
      double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
      double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2]));
      
      if(angle<_coneAngle)
	{
	  if(angle>OpAngleMax)
	    OpAngleMax=angle;
	  tau.push_back(Qvec[s]);
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    streamlog_out(DEBUG) << "Adding Q: "<<track->getType()<<"\t"<<track->getEnergy()<<"\t"<<p<<"\t"<<theta<<"\t"<<phi<<std::endl;
	  Etau+=Qvec[s]->getEnergy();
	  //combine to new momentum
	  for(int i=0;i<3;i++){
	    pvec_tau[i]=pvec_tau[i]+track->getMomentum()[i];
	  }
	  Qvec.erase(iterQ);
	  s--;
	}
      else
	iterQ++;
    }
  //assign neutral particles to tau candidate
  std::vector<ReconstructedParticle*>::iterator iterN=Nvec.begin();
  for (unsigned int s=0; s<Nvec.size() ; s++ )
    {
      ReconstructedParticle *track=Nvec[s];
 
      double *pvec=(double*)track->getMomentum();
      double angle=acos((pvec[0]*pvec_tau[0]+pvec[1]*pvec_tau[1]+pvec[2]*pvec_tau[2])/
			(sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2])*
			 sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]+pvec_tau[2]*pvec_tau[2])));
      double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
      double p=sqrt(pt*pt+pvec[2]*pvec[2]);
      double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
      double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2]));
      
      if(angle<_coneAngle)
	{
	  if(angle>OpAngleMax)
	    OpAngleMax=angle;
	  tau.push_back(Nvec[s]);
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
            streamlog_out(DEBUG) << "Adding N: "<<track->getType()
                                 <<"\t"<<track->getEnergy()
                                 <<"\t"<<p
                                 <<"\t"<<theta
                                 <<"\t"<<phi<<std::endl;

	  Etau+=Nvec[s]->getEnergy();
	  //combine to new momentum
	  for(int i=0;i<3;i++){
	    pvec_tau[i]=pvec_tau[i]+track->getMomentum()[i];
	  }
	  s--;
	  Nvec.erase(iterN);
	}
      else 
	iterN++;
    }
  tauvec.push_back(tau);
  return false;
}

void TauFinder::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TauFinder::end(){ 
  anatree->Write();
  rootfile->Write();
  delete anatree;
  delete rootfile;
 

}
