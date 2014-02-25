#include <TMath.h>

#include "TSharc.h"


ClassImp(TSharc)

int TSharc::fCfdBuildDiff = 5;


TSharc::TSharc()	{	}

TSharc::~TSharc()	{	}





void	TSharc::BuildHits(Option_t *opt)	{

  //  after the data has been taken from the fragement tree, the data
  //  is stored/correlated breifly in by the tsharcdata class - these 
  //  function takes the data out of tsharcdata, and puts it into the 
  //  the tsharchits class.  These tsharchit objects are how we access
  //  the data stored in the tsharc branch in the analysis tree. 
  //
  //  pcb.
  //

  TSharcHit sharchit;
  //int fCfdBuildDiff = 5; // largest acceptable time difference between events (clock ticks)  (50 ns)

	std::map<std::tuple<int,int,int>,int> hitmap;

  for(int i=0;i<GetMultiplicityFront();i++)	{	
  	
    for(int j=0;j<GetMultiplicityBack();j++)	{	

      if(GetFront_DetectorNbr(i) == GetBack_DetectorNbr(j))	{ //check if same detector
		
				if(abs(GetFront_TimeCFD(i)-GetBack_TimeCFD(j)) > fCfdBuildDiff)
					continue; // ensure there s front-back time correlation to protect against noise
	
				if(abs((float)GetFront_Charge(i)-(float)GetBack_Charge(j))/((float)GetFront_Charge(i) + (float)GetBack_Charge(j)) > 0.201)
					continue; // ensure the difference between the detectred front charge
							  // and back charge are with in 20% of each other.
		
				sharchit.SetDetector(GetFront_DetectorNbr(i));
	
				sharchit.SetDeltaE(GetFront_Energy(i)); 
				sharchit.SetFrontCharge(GetFront_Charge(i));
				sharchit.SetBackCharge(GetBack_Charge(j));
				sharchit.SetDeltaT(GetFront_Time(i)) ;  		//cheak time allignment;

				sharchit.SetFrontCFD(GetFront_TimeCFD(i));
				sharchit.SetBackCFD(GetBack_TimeCFD(j));			

				sharchit.SetPixel(GetFront_StripNbr(i),GetBack_StripNbr(j));
				//sharchit.SetDeltaCfd(GetFront_TimeCFD(i));		//cheak time allignment;

				sharchit.SetPosition(TSharc::GetPosition(GetFront_DetectorNbr(i),GetFront_StripNbr(i),GetBack_StripNbr(j)));
																				
				sharc_hits.push_back(sharchit);																				

				hitmap.insert({std::make_tuple((int)GetFront_DetectorNbr(i),(int)GetFront_StripNbr(i),(int)GetBack_StripNbr(j)),sharc_hits.size()-1});


      }
    }
  }

   if(opt == "CLEAN")	{
	  //printf("\nINSIDE CLEAN IF!!!!!!!!\n");
	  //fflush(stdout);
  	  std::set<int> used_hits;	
      std::map<std::tuple<int,int,int>,int>::iterator iter1;	
      std::map<std::tuple<int,int,int>,int>::iterator iter2;	

      for(iter1 = hitmap.begin(); iter1 != hitmap.end(); iter1++)	{
         for(iter2 = std::next(iter1); iter2 != hitmap.end(); iter2++)	{
			printf("\titer1 < %02i | %02i | %02i >\n",std::get<0>(iter1->first),std::get<1>(iter1->first),std::get<2>(iter1->first));
			printf("\titer2 < %02i | %02i | %02i >\n",std::get<0>(iter2->first),std::get<1>(iter2->first),std::get<2>(iter2->first));
			printf("\t=============================\n");
			if(std::get<0>(iter1->first) != std::get<0>(iter2->first))  ///  check detector numbers;
         		break;
         	if( abs(std::get<1>(iter1->first) - std::get<1>(iter2->first))<2) {  ///adjacent front strips.
				if( (used_hits.count(iter1->second) == 0) &&  (used_hits.count(iter2->second) == 0) ) {
					used_hits.insert(CombineHits(&sharc_hits.at(iter1->second),&sharc_hits.at(iter2->second),iter1->second,iter2->second));
         		}
			}
         	else if( abs(std::get<2>(iter1->first) - std::get<2>(iter2->first))<2) {  ///adjacent backstrips strips.
				if( (used_hits.count(iter1->second) == 0) &&  (used_hits.count(iter2->second) == 0) ) {
					used_hits.insert(CombineHits(&sharc_hits.at(iter1->second),&sharc_hits.at(iter2->second),iter1->second,iter2->second));
         		}
			}
			else	{
				break;
			}
         }
      }
      RemoveHits(&sharc_hits,&used_hits);
  }




  for(int k=0;k<GetMultiplicityPAD();k++)	{	
	//int oldtotalhits = TSharc::totalhits;
    for(int l=0;l<sharc_hits.size();l++)	{
      if(GetPAD_DetectorNbr(k) == sharc_hits.at(l).GetDetectorNumber())	{ //check if same detector
				
				sharc_hits.at(l).SetPadE(GetPAD_Energy(k));
				sharc_hits.at(l).SetPadCharge(GetPAD_Charge(k));
				sharc_hits.at(l).SetPadT(GetPAD_Time(k));
				//sharc_hits.at(l).SetPadCFD(GetPAD_TimeCFD(k));
      }
    }
  }
}


//TVector3 TSharc::GetPosition(TSharcHit *hit)	{
//
//   The vector3 must be stored in the TSharcHit class to 
//   allow for weighted avgs of differenet pixel to be calculated
//   and stored.
//
//}


int TSharc::CombineHits(TSharcHit *hit1,TSharcHit *hit2,int position1,int position2 )	{
	// used in the build hits routine to combine to hits 
	// in adcent pixels into one hit.

	if(!hit1 || !hit2)
		return 0xffffffff;

	float hit1_weight = hit1->GetFrontCharge()/(hit1->GetFrontCharge() + hit2->GetFrontCharge());
	if(hit1_weight < 0.051)
		return position1;
	float hit2_weight = hit2->GetFrontCharge()/(hit1->GetFrontCharge() + hit2->GetFrontCharge());
	if(hit2_weight < 0.051)
		return position2;

	if(hit1_weight > hit2_weight)	{
		hit1->SetFrontCharge(hit1->GetFrontCharge() + hit2->GetFrontCharge());		
		hit1->SetBackCharge(hit1->GetBackCharge() + hit2->GetBackCharge());		
		TVector3 newposition;
		newposition.SetX(hit1->GetPosition().X()*hit1_weight + hit2->GetPosition().X()*hit2_weight);
		newposition.SetY(hit1->GetPosition().X()*hit1_weight + hit2->GetPosition().X()*hit2_weight);
		newposition.SetZ(hit1->GetPosition().X()*hit1_weight + hit2->GetPosition().X()*hit2_weight);
		hit1->SetPosition(newposition);
		return position2;
	}
	else {
		hit2->SetFrontCharge(hit1->GetFrontCharge() + hit2->GetFrontCharge());		
		hit2->SetBackCharge(hit1->GetBackCharge() + hit2->GetBackCharge());		
		TVector3 newposition;
		newposition.SetX(hit1->GetPosition().X()*hit1_weight + hit2->GetPosition().X()*hit2_weight);
		newposition.SetY(hit1->GetPosition().X()*hit1_weight + hit2->GetPosition().X()*hit2_weight);
		newposition.SetZ(hit1->GetPosition().X()*hit1_weight + hit2->GetPosition().X()*hit2_weight);
		hit2->SetPosition(newposition);
		return position1;
	}
}

void TSharc::RemoveHits(std::vector<TSharcHit> *hits,std::set<int> *to_remove)	{

	std::set<int>::reverse_iterator iter;
//	for(iter = to_remove->begin(); iter != to_remove->end(); iter++)	{
//		if(*iter == 0xffffffff)
//			continue;
//		hits->erase(std::remove(hits->begin(),hits->end(),*iter ),hits->end());
//	}
	for(iter= to_remove->rbegin(); iter != to_remove->rend(); iter++)	{
		if(*iter == 0xffffffff)
			continue;
		hits->erase(hits->begin()+*iter);

	}
}










void TSharc::Clear(Option_t *option)	{

  //cout << "clearing " << endl;
  ClearData();
  sharc_hits.clear();
  //cout <<" size: " << sharc_hits.size() << endl;
  return;
}

void TSharc::Print(Option_t *option)	{
  printf("not yet written...\n");
  return;
}



//double TSharc::GetDeadLayer(TSharcHit *hit)	{
//	
//
//}




TVector3 TSharc::GetPosition(int detector, int frontstrip, int backstrip)	{
  int FrontDet = detector;
  int FrontStr = frontstrip;
  int BackDet  = detector;
  int BackStr  = backstrip;


  TVector3 position;
  double stripBpitch = 48.0/48;
  double stripFpitch = 72.0/24;
  double Zoffset = 7.;
  double Yoffset = 72.0/2;
  double x = 0;
  double y = 0;
  double z = 1;
  if(FrontDet>4 && FrontDet<9){ //forward box
    z = stripBpitch*(BackStr+0.5)+Zoffset;
    y = Yoffset;
    x = stripFpitch*(FrontStr-12+0.5);
    position.SetXYZ(x,y,z);
  }
  else if(FrontDet>8 && FrontDet<13){ //backward box
    z = stripBpitch*(BackStr+0.5)+Zoffset;
    z= -z;
    y = Yoffset;
    x = stripFpitch*(FrontStr-12+0.5);
    position.SetXYZ(x,y,z);
  }
  else{
    double outerradius =41;
    double ringpitch = (41.-9.)/16;
    double segmentpitch = 81.6/24;
    double z = 66.0;
    double rho = outerradius - (FrontStr+0.5)*ringpitch;
    double phi = (BackStr+0.5)*segmentpitch*TMath::Pi()/180.;
    phi+=(4-FrontDet)*TMath::Pi()/2.;
    position.SetXYZ(1.,1.,1.);
    position.SetPerp(rho);
    position.SetPhi(phi);
    if(FrontDet<5)
      position.SetZ(z);
    else
      position.SetZ(-z);
  }
  return position;
}

double TSharc::GetDeadLayer(int detector, int frontstrip, int backstrip)		{
    TVector3 position = TSharc::GetPosition(detector,frontstrip,backstrip);

    double deadlayerthicknessinput[16] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};  // note: make a .sharc_rc to put stuff like this.

    if(detector>4 && detector<13 ){ // BOX
        return deadlayerthicknessinput[detector-1]/(  TMath::Sin(position.Theta())*TMath::Cos(position.Phi() - (TMath::Pi()/2)) );
	}
    else if( detector <= 4 || detector >= 13 ){ // QQQ
        return deadlayerthicknessinput[detector-1]/(TMath::Cos(position.Theta()));
    }

}










