#include <iostream>
#include <list>
#include <math.h>
#include <stdlib.h>
#include "stems.h"
#include "nr_RSG.h"


#define NEW_BALANCE // Flag to turn on random placement of new stems in the list (enforces detailed balance for all microstates)

//===About 9kBT at 20 segs( for ANG_TOL = 4PI)=====================
//==Nuc time = 3.138215e+04 tau ===================================
#define EBULK 0.0//1.9
#define ESURFACE 0.0//1.9
#define AR (1.33 * 1.0)

//===About 8kBT at 30 segs( for ANG_TOL = 2PI)=====================
//Can compensate for change in angular tollerance by adding -ln( change in angTol)
//on to EBULK. This preserves the free energy landscape but changes the kinetics!
//#define EBULK 0.03
//#define ESURFACE 0.27//0.22 0.25



#define PI 3.141592654
#define ANG_TOL (4 * PI)
//#define ANG_TOL (2 * PI)






#define TOTAL_SEGS 500 //no more than 270
#define TOTAL_STEMS 500
//#define LONGEST_STEM 2
#define MIN_STEMS 1
#define MIN_SEGS 1
//#define NO_TOP 1// don't allow top or bottom addition or removal!
//#define NO_MOVE 1 // dont make any moves!!!

#define FINISH_TIME 1.0e5

static void firstPoint( list < StemPtr >  *cry);
static void newPoint(list < StemPtr > *cry);
static void addStem(list < StemPtr > *cry, int addSpecies);

extern void printAll( list < StemPtr > cry);
void printInt( int x);
int countSegments( list < StemPtr > crystal );
void catiousTesting( list < StemPtr > cry );


void calculateDegeneracy( double degen[TOTAL_SEGS+1][TOTAL_STEMS+1] );

double surfaceAreaEllipse( int NT, int Ns);
double freeEnergy( int NT, int Ns);
static double computeAspectRatio( double Nt, double Ns);



static double sumRates( list < StemPtr > cry);
static int pickAndMakeMove( list < StemPtr > *cry, double rateSum );
//static void pickMove( list < StemPtr > cry, double rateSum, 
//		      int *event,list<StemPtr>::iterator *evItr);
//static void makeMove( list < StemPtr > *cry, int move, 
//		      list<StemPtr>::iterator *evItr);

static double sideAreaFunction( int Ns);
double degeneracyFactor( int NT, int Ns );

double getDegPoint( double degen[TOTAL_SEGS+1][TOTAL_STEMS+1], int seg, int stem );

double facRatio( long a, long b);

double factorial( long N);
double getDegen(int T, int S);


void zapCrystal(list < StemPtr > *cry);





void plotLandscape();
void findNucTime();
void findNucFrac();


double computeTotalFreeEnergy( int Nt, double degen[TOTAL_SEGS+1][TOTAL_STEMS+1]);
template <class T> inline const T& myMin ( const T& a, const T& b );
template <class T> inline const T& myMax ( const T& a, const T& b );


int nStems, nSegs;
double eBulk, eSurface;
double addStemRate_global; //needed to save recomputing this value

long seed;



int main()
{
  int i,j, counter=0;
  int max=9;
  

 

  plotLandscape();
  //findNucTime ();
  //findNucFrac();

  return 0;
}


  
void findNucTime()
{

  //*******CHECK THE SIDE AREA FUNCTION IS WHAT YOU EXPECT IT TO BE!!!!!******************************************

#ifdef TOTAL_SEGS
   list<StemPtr> crystal;
   long i;
   int j,event, runNumber;
   double rateSum, sTime, timeIncre,  sum;
   list<StemPtr>::iterator iL, eventItr;
   double accTime=0.0, aspRatio_acc=0.0;
   char name[100];
   char mastername[100];
   FILE *cfPtr;

   printf("Enter the Filename :");
   scanf("%s",name);

   sprintf(mastername,"Data/%s.dat",name);
   
   

   //=========initialise variables============================
   //eBulk = 0.1;   eSurface = 0.20; // gives an energy barrier of about 2.3 at 12 segments
   eBulk = EBULK;   eSurface = ESURFACE; 

   //eBulk = 0.0;eSurface = 0.0;
   

   seed = -time(NULL);
   

   for( runNumber =1; runNumber <= 1000000; runNumber++){
   
     //==Initialise crystal==========================================
     firstPoint( &crystal );
     //addStem( &crystal , 1 );

     
     sTime = 0.0;
     catiousTesting( crystal );

     
     while (nSegs < TOTAL_SEGS){
       
     //*****************************************************
     //THESE TWO FUNCTIONS SHOULD AWAYS BE CALLED TOGETHER
       rateSum = sumRates( crystal );
       event = pickAndMakeMove( &crystal, rateSum);
       //******************************************************
       
       
       //==========Sort out time and check crystal===============
       timeIncre = -log(ran1( &seed))/ rateSum ;
       sTime += timeIncre;
       catiousTesting( crystal );
     }

     accTime += sTime;
     aspRatio_acc += computeAspectRatio( nSegs, nStems);
     
     //====Write data and close file================
     //(open/writing/closing the file appears to have a negligible effect on the run time)
     cfPtr = fopen (mastername,"a");
     if (cfPtr==NULL){
       printf("Can't open file %s\n",mastername);
       abort();
     }
     fprintf(cfPtr,"%e\n",sTime);
     fclose(cfPtr);
     

     if(fmod(runNumber,100)==0)
       printf("%d %f %e %f %d\n",runNumber, 1.0*TOTAL_STEMS,accTime/(1.0*runNumber)
	      , aspRatio_acc/(1.0*runNumber), runNumber);
     

     
     //delete the crystal from memory.
     zapCrystal( &crystal );
   }

#endif

   return;
}

void findNucFrac()
{

  //*******CHECK THE SIDE AREA FUNCTION IS WHAT YOU EXPECT IT TO BE!!!!!**************************************************************


#ifdef TOTAL_SEGS
   list<StemPtr> crystal;
   long i;
   int j,event, runNumber;
   double rateSum, sTime, timeIncre,  sum;
   list<StemPtr>::iterator iL, eventItr;
   int accEvents=0;
   bool nuc;

   //FILE *cfPtr;
   

   //=========initialise variables============================
   //eBulk = 0.1;   eSurface = 0.20; // gives an energy barrier of about 2.3 at 12 segments
   eBulk = EBULK;   eSurface = ESURFACE; 

   //eBulk = 0.0;eSurface = 0.0;
   

   seed = -time(NULL);
   

   for( runNumber =1; runNumber <= 1000000000; runNumber++){
   
     //==Initialise crystal==========================================
     firstPoint( &crystal );
     

     
     sTime = 0.0;
     catiousTesting( crystal );

     nuc = false;
     while (nSegs < TOTAL_SEGS && sTime<FINISH_TIME){
       
       

     //*****************************************************
     //THESE TWO FUNCTIONS SHOULD AWAYS BE CALLED TOGETHER
       rateSum = sumRates( crystal );
       event = pickAndMakeMove( &crystal, rateSum);
       //******************************************************
       
       
       //==========Sort out time and check crystal===============
       timeIncre = -log(ran1( &seed))/ rateSum ;
       sTime += timeIncre;
       catiousTesting( crystal );

       if( nSegs == TOTAL_SEGS && sTime < FINISH_TIME ) nuc=true;

     }

     
     if( nuc ) accEvents ++;

     if(fmod(runNumber,1000)==0)
       printf("%f %e %d\n",  FINISH_TIME, 1.0*accEvents/(1.0*runNumber), runNumber);

     

     
     //delete the crystal from memory.
     zapCrystal( &crystal );
   }

#endif

   return;
}


  
void plotLandscape()
{
#ifdef TOTAL_SEGS
   list<StemPtr> crystal;
   long i,  segsLastStep, stemsLastStep;
   int j,k,event;
   double rateSum, sTime, timeIncre,  sum;
   //   double F[TOTAL_SEGS+1]={0.0}, partFn;
   list<StemPtr>::iterator iL, eventItr;
   double accTime[TOTAL_SEGS+1]={0.0};

   double accTime2[TOTAL_SEGS+1][TOTAL_STEMS+1]={0.0};
   double degen[TOTAL_SEGS+1][TOTAL_STEMS+1]={0.0};
   long int state211=0, state121=0, state112=0;
   long int state22=0, state13=0, state31=0;



   FILE *cfPtr;
   

   

   
   //=========initialise variables============================
   eBulk = EBULK;   eSurface = ESURFACE; // gives an energy barrier of about 2.3 at 12 segments
   //eBulk = 0.0;
   //eSurface = 0.0;
   

   seed = -time(NULL);
   

   
   calculateDegeneracy( degen );
  

   cfPtr = fopen("Calc.dat","w");

   for(i=1; i<=TOTAL_SEGS ; i++){
     fprintf(cfPtr,"%f %e\n",1.0*i, computeTotalFreeEnergy( i , degen ) );
   }

   fclose( cfPtr);

   //return;
   

   //==Initialise crystal==========================================
   firstPoint( &crystal );

   //(*(crystal.begin()))->top ++;
   //   addStem( &crystal , 1 );
   printAll( crystal );
   rateSum = sumRates( crystal );
   printf("rateSum=%f\n",rateSum);
   //return;
   
   //addStem( &crystal , 1 );
   //segsLastStep = 1;

  



   sTime = 0.0;
   catiousTesting( crystal );


   for(i=1; i<= 1000000000;i++){

     //printAll( crystal );

     segsLastStep = nSegs; //must log previous step to get free energy correct
     stemsLastStep = nStems;
     
     
     

     //*****************************************************
     //THESE TWO FUNCTIONS SHOULD AWAYS BE CALLED TOGETHER
     rateSum = sumRates( crystal );
     event = pickAndMakeMove( &crystal, rateSum);
     //******************************************************
     

     //==========Sort out time and check crystal===============
     timeIncre = -log(ran1( &seed))/ rateSum ;
     sTime += timeIncre;
     catiousTesting( crystal );


     


     //===========Log time in state before timestep===============
     accTime[segsLastStep] += timeIncre;
     accTime2[segsLastStep][stemsLastStep] += timeIncre;


     
     if( nSegs==4 && nStems==3){
       iL= crystal.begin();
       //printf("%d\n",(*iL)->top);

       if( (*iL)->top - (*iL)->bot==1) {
	 state211++;
	 //printf("State 12\n");
	 //printAll(crystal);
	 //printf("\n");
       }else{
	 ++iL;
	 if( (*iL)->top - (*iL)->bot==1) {
	   //printf("State 21\n");
	   //printAll(crystal);
	   //printf("\n");
	   state121++;
	 }else{
	   state112++;
	 }
       }
     }

     if( nSegs==4 && nStems==2){
       iL= crystal.begin();
       if( (*iL)->top - (*iL)->bot==0) state13++;
       if( (*iL)->top - (*iL)->bot==1) state22++;      
       if( (*iL)->top - (*iL)->bot==2) state31++;
     }

     
     /*
     printf("%d %d\n",nStems, nSegs);
     printAll( crystal );
     */

     
     //==========Output data======================================
     if( fmod(i,100000) == 0){
       
       sum = 0.0;
       for(j=1; j<=TOTAL_SEGS ; j++)	 sum += accTime[j];


       cfPtr = fopen("data.dat","w");

       for( j=1; j<=TOTAL_SEGS; j++){
	 //printf("%d) %f %f\n",j, accTime[j]/sTime, F[j]/partFn);
	 if( accTime[j]>0){
	   fprintf(cfPtr,"%f %f %f\n",1.0*j, -log(accTime[j]/sTime)+log(accTime[1]/sTime),
		   accTime[j]/sTime);
	 }
       }

   
       fclose( cfPtr);
       printf("%ld %f %f %f %f\n",i, 1.0*state211/(1.0*state121) , 1.0*state112/(1.0*state121), 
	      1.0*state13/(1.0*state31), 1.0*state22/(1.0*state31));

     }
     

   }
   
   printAll( crystal );

#endif

   return;
}


static double sumRates( list < StemPtr > cry)
{
  //determines rate of moves and sums over all moves
  list<StemPtr>::iterator iL;
  double sum, currentFE, dF_addTop, dF_removeTop, dF_addStem, dF_removeStem,
    removeStemScale;


  


  //determine changes in free energy
  currentFE = freeEnergy( nSegs, nStems);
  dF_addTop = exp(-freeEnergy( nSegs +1, nStems) + currentFE);
  dF_removeTop = exp(-freeEnergy( nSegs -1, nStems) + currentFE);
  dF_addStem = exp(-freeEnergy( nSegs +1, nStems+1) + currentFE);
  dF_removeStem = exp(-freeEnergy( nSegs -1, nStems-1) + currentFE);

  //make new value of addStem available globally
  addStemRate_global = sideAreaFunction(nStems) * myMin( 1.0, dF_addStem);

  //compute scale factor for stem removal
  removeStemScale = sideAreaFunction( nStems -1 ) / (1.0*nStems) * degeneracyFactor(nSegs, nStems);
#ifdef MIN_STEMS
  if(nStems == MIN_STEMS) removeStemScale = 0.0;
#endif


  //printf("Current = %f, plus stem = %f, add seg=%f\n",currentFE, 
  //	 freeEnergy( nSegs +1, nStems+1), freeEnergy(nSegs+1,nStems));


  sum = addStemRate_global * ANG_TOL/(4 * PI); // there are sqrt(nStems) possible side addition moves

  //if it's the only segment then prevent removal
  if(nSegs == 1 ){
    iL=cry.begin();
    (*iL)->addTop = myMin( 1.0, dF_addTop) * ANG_TOL/(4 * PI);
    (*iL)->addBot = myMin( 1.0, dF_addTop) * ANG_TOL/(4 * PI);
    (*iL)->removeBot = 0.0;
    (*iL)->removeTop = 0.0;
    
    sum += (*iL)->addTop + (*iL)->addBot + (*iL)->removeTop + (*iL)->removeBot;

  }else{
        
    for(iL=cry.begin(); iL != cry.end(); ++iL){
      (*iL)->addTop = myMin( 1.0, dF_addTop) * ANG_TOL/(4 * PI);
      (*iL)->addBot = myMin( 1.0, dF_addTop) * ANG_TOL/(4 * PI);
      
      // if the top equals the bottom then just count the removal move once
      if( (*iL)->top == (*iL)->bot ) {
	(*iL)->removeTop = myMin(1.0,  dF_removeStem) * removeStemScale;
	(*iL)->removeBot = 0.0;
      }else{
	(*iL)->removeTop = myMin( 1.0, dF_removeTop) * degeneracyFactor(nSegs, nStems);;
	(*iL)->removeBot = myMin( 1.0, dF_removeTop) * degeneracyFactor(nSegs, nStems);;
      }
      
      sum += (*iL)->addTop + (*iL)->addBot + (*iL)->removeTop + (*iL)->removeBot;
    }
  }

  return sum;
}

static int pickAndMakeMove( list < StemPtr > *cry, double rateSum )
{

  //Picks the event according to rate probabilities and then returns the event number, along
  //with an pointer to the stem the event occurs on.
  
  double ranNum, rateAcc=0;
  int event;
  list<StemPtr>::iterator iP;
  
  
  //fix random number
  ranNum = rateSum * ran1(&seed);
  
  
  //cycle through all possible events
  event = 1;
  rateAcc += addStemRate_global * ANG_TOL/(4 * PI);
  if( ranNum < rateAcc ){
#ifdef TOTAL_SEGS
    if( nSegs == TOTAL_SEGS) return event;// care!!
#endif
#ifdef NO_MOVE
    return event;
#endif
    addStem( cry , 1 ); // add stem also incremens nSegs etc
    //printf("Event 1\n");
    return event;
  }
  
  
  for(iP=(*cry).begin(); iP != (*cry).end(); ++iP){
    
    //==================add top=============================
    (event) ++;
    rateAcc += (*iP)->addTop;
    if( ranNum < rateAcc ){
#ifdef TOTAL_SEGS
      if( nSegs == TOTAL_SEGS) return event;// care!!
#endif

#ifdef NO_TOP
      return event;
#endif

#ifdef LONGEST_STEM
      if( (*iP)->top-(*iP)->bot +1 == LONGEST_STEM) return event;
#endif


      (*iP)->top++;
      nSegs ++;
      //remember to reassign Z!!!!!
      //printf("pickMove %d %d %d %d\n",event, (*iP)->top, (*iP)->bot, 
      //(*iP)->species);
      
      return event;
    }
    
    
    //========================remove top==========================
    (event) ++;
    rateAcc += (*iP)->removeTop;
    if( ranNum < rateAcc ){
#ifdef NO_MOVE
      return event;
#endif

#ifdef MIN_SEGS
    if( nSegs == MIN_SEGS) return event;// care!!
#endif

      if( (*iP)->top == (*iP)->bot ){
	free( *iP );
	iP = (*cry).erase( iP );
	nStems --;
      }else{
#ifdef NO_TOP
      return event;
#endif
      (*iP)->top --;
      }
      nSegs --;
      
      //printf("pickMove rt %d %d %d %d\n",event, (*iP)->top, (*iP)->bot,    (*iP)->species);
      return event;
    }
    
    
    //====================add Bot=====================================
    (event) ++;
    rateAcc += (*iP)->addBot;
    if( ranNum < rateAcc ){
#ifdef TOTAL_SEGS
      if( nSegs == TOTAL_SEGS) return event;// care!!
#endif

#ifdef NO_TOP
      return event;
#endif


#ifdef LONGEST_STEM
      if( (*iP)->top-(*iP)->bot +1 == LONGEST_STEM) return event;
#endif

      //printf("pickMove %d %d %d %d\n",event, (*iP)->top, (*iP)->bot,  (*iP)->species);
      (*iP)->bot --;
      nSegs ++;
      return event;
    }
    
    //================================remove Bot=========================
    (event) ++;
    rateAcc += (*iP)->removeBot;
    if( ranNum < rateAcc ){
      
      if( ranNum < rateAcc ){

#ifdef NO_TOP
      return event;
#endif

#ifdef MIN_SEGS
    if( nSegs == MIN_SEGS) return event;// care!!
#endif


	if( (*iP)->top == (*iP)->bot ){
	  printf("This should never happen!\n");
	  abort();
	}else{
	  (*iP)->bot ++;
	}
	nSegs --;
	//printf("pickMove %d %d %d %d\n",event, (*iP)->top, (*iP)->bot,      (*iP)->species);
	return event;
      }
      
      
    }
  }
  
  //never reach here
  printf("Error\n");
  abort();
  
  
  return -1;
}

/*
static void makeMove( list < StemPtr > *cry, int move, 
		      list<StemPtr>::iterator *evItr)
{
  //Receives an integer event number (move) and an iterator pointing
  // to the stem of the event.
  //Carrys out the select evented on the indicated stem (pointer!)

  int stemNumber, moveNumber;
  //list<StemPtr>::iterator iL,;


  printf("Making move %d\n",move);
  printf("Top =%d bot=%d\n",(*(*evItr))->top, (*(*evItr))->bot);


  if( move ==1 ){
    printf("Move 1 caught\n");
    addStem( cry , 1 );
    return;
  }

  stemNumber =  (int)ceil((move-1) /4.0);// regular integer division rounds down so add 1 to round up
  moveNumber = move - 1 - (stemNumber-1)*4;

  printf("move = %d, stem=%d, mNumber=%d\n",move, stemNumber, moveNumber);
  printf("Top =%d bot=%d\n",(*(*evItr))->top, (*(*evItr))->bot);
  
  if( moveNumber == 1){
    //#ifdef    
    (*(*evItr))->top ++;
    return;
  }


  return;
}
*/


extern void firstPoint(list < StemPtr > *cry)
{
  StemPtr newPtr;

  newPtr = new stem;
 
   (*cry).push_front( newPtr );

   newPtr-> top =50;
   newPtr-> bot =50;
   newPtr-> species =1;

   newPtr-> Ztop = 0;
   newPtr-> Zbot = 0;


   nStems = 1;
   nSegs =1;



  return;
}


static void addStem(list < StemPtr > *cry, int addSpecies)
{
  int i;
  StemPtr newPtr;
  list<StemPtr>::iterator it;
  int newPosition=0;



#ifdef TOTAL_STEMS
  if( nStems == TOTAL_STEMS) return; //CARE!!!!!!!
#endif



#ifdef NEW_BALANCE
  //generate random number between 0 and nStems
  newPosition = int(ran1( &seed) * (nStems+1) );

  if( newPosition> nStems+1 || newPosition<0){
      printf("I have failed\n");
      abort();
    }
#endif



  //newPosition = myMin(2,nStems);
  //printf("New pos is after point %d\n",newPosition);

  newPtr = new stem;
  
  if( newPosition == 0){
    (*cry).push_front( newPtr );

  }else{
    if( newPosition == nSegs){
      (*cry).push_back( newPtr );

    }else{

      it = (*cry).begin();
      for(i=1; i<=newPosition ; i++) ++it;
      (*cry).insert (it , newPtr);

    }
  }

  
   newPtr-> top =50;
   newPtr-> bot =50;
   newPtr-> species = addSpecies;


   
   newPtr-> Ztop = 0;
   newPtr-> Zbot = 0;


   nStems ++;
   nSegs ++;

   return;
}




double freeEnergy( int NT, int Ns)
{

  if( NT == 0 ) return 0;

  return  -eBulk*NT  + eSurface * surfaceAreaEllipse( NT, Ns ); 
  //return 1;
}

double surfaceAreaEllipse( int NT, int Ns)
{  
  double aspRatioSq, ellipticity;

  aspRatioSq =  AR * AR * NT *NT / Ns/Ns/Ns;

  //printf("%d %f\n",Ns,aspRatioSq);

  if( fabs(aspRatioSq -1.0) <1e-5 ) return 4.0 * Ns;

  if( aspRatioSq >1.0 ){

    ellipticity = sqrt( 1.0 - 1.0/aspRatioSq );
    return 2.0 * Ns   +   2.0* AR *NT / sqrt(1.0*Ns) / ellipticity * asin( ellipticity ) ;

  }else{

    ellipticity = sqrt( 1.0- aspRatioSq);    
    //ellipticity = acos( length/width);
    return  2.0*Ns   +   AR * AR * NT*NT/(1.0*Ns*Ns) / ellipticity * log( (1.0 + ellipticity)/ (1.0-ellipticity) );

  }
 
  abort();
  return 1;

}

void catiousTesting( list < StemPtr > cry )
{

  if( countSegments( cry) != nSegs ){
    printAll( cry);
    printf("%d !=  %d\n",  countSegments( cry) , nSegs);
    abort();
  }  
  return;
}

int countSegments( list < StemPtr > cry )
{  
  list<StemPtr>::iterator iL;
  int sum = 0;

  for(iL=cry.begin(); iL != cry.end(); ++iL){
    sum += (*iL)->top - (*iL)->bot +1;
  }

  return sum;
}



extern void printAll( list < StemPtr > cry)
{
  int i;
  list<StemPtr>::iterator iL;
  double rateSum;

  rateSum = sumRates( cry );

  printf("nSegs = %d rateSum=%f\n",nSegs, rateSum);

  printf("        sp  top bot Zt  Zb  T+       T-       B+       B-\n");

  i=1;
  for(iL=cry.begin(); iL != cry.end(); ++iL){

    printInt(i);
    printf(") ");
    printInt( (*iL)-> species );
    printInt( (*iL)-> top );
    printInt( (*iL)-> bot );
    printInt( (*iL)-> Ztop );
    printInt( (*iL)-> Zbot );
    
	     

    printf("  %.2e %.2e %.2e %.2e\n", (*iL)->addTop, (*iL)->removeTop, (*iL)->addBot, (*iL)->removeBot );
    i++;
  }

  printf("\n");

  return;
}
   
void printInt( int x)
{
  if( x< 10 ){
    printf("  %d ",x);
    return;
  }

  if( x< 100 ){
    printf(" %d ",x);
    return;
  }

  printf("%d ",x);
  return;
}


template <class T> inline const T& myMin ( const T& a, const T& b ) {
  return (a<b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}

template <class T> inline const T& myMax ( const T& a, const T& b ) {
  return (a>b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}


static double sideAreaFunction( int Ns)
{
  //return 0.0;
  //return 1.0*Ns*Ns;
  return 3.545 * sqrt( 1.0* Ns); //3.545 = 2 * sqrt(pi)
  //return 1.0;
}

double degeneracyFactor( int NT, int Ns )
{
  //a testing function that lets me play around with different ways to cancel out degeneracy as 
  //the number of segments increases

  

  if( NT==-3)
    return 1.0*Ns;
  else
    return 1.0;
}

void calculateDegeneracy( double degen[TOTAL_SEGS+1][TOTAL_STEMS+1] )
{
  int row, col;
  int step;
  


  for(step=1; step<=TOTAL_SEGS ; step++){
    
    //seed this row
    degen[step][1] = 1;

    //move along bottom row
    for(col=2; col<=step; col++) 
      degen[step][col] = degen[step-1][col] + degen[step-1][col-1];
    
    //seed final corner (all other points in this column are zero
    degen[step][step] = 1;
      

    
  }


  return;
}

  
double factorial( long N)
{
  double i, prod=1;
 

  for(i=1; i<=N ; i++)
    prod = prod *i;
  
  
  return prod;
}

double facRatio( long a, long b)
{
  //returns a!/b! where a >= b

  double prod =1;
  int i;

  if(a==b) return 1;

  if(a==b+1) return a;

  for(i=b+1; i<=a; i++)
    prod = prod * i;

  return prod;
}
  


double getDegen(int T, int S)
{

  if (S-1> T-S) 
    return facRatio(T-1, S-1)/ factorial(T-S);
  else
    return facRatio(T-1, T-S) / factorial(S-1);
}

double computeTotalFreeEnergy( int Nt, double degen[TOTAL_SEGS+1][TOTAL_STEMS+1])
{
  double Zn=0.0;
  int stem;
  
  for(stem=1; stem<=Nt ; stem++){
    Zn += degen[Nt][stem] * exp( - freeEnergy(Nt,stem) );
    //printf("w=%e B=%e FE=%e\n",degen[Nt][stem] ,exp( - freeEnergy(Nt,stem) ) , freeEnergy(Nt,stem) );

  }

  return - log( Zn ) - freeEnergy(1,1);
}


void zapCrystal(list < StemPtr > *cry)
{
  list<StemPtr>::iterator iL;

  //free up the individual structures
  for(iL=(*cry).begin(); iL != (*cry).end(); ++iL)
    free( *iL );
  
  //clear the list of pointers
  (*cry).clear();

  //reset stem and segment counters
  nSegs = 0;
  nStems = 0;

  return;
}
  
static double computeAspectRatio( double Nt, double Ns)
{
  return AR * Nt/ Ns/ sqrt(Ns); 
}  
