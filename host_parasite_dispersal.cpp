/*
	Copyright (C) 2020  Jhelam N. Deshpande
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

//============================================================================


#include <cstdlib>					// standard C library
#include <iostream>					// standard input/output functions library
#include <fstream>					// file stream library
#include <string>					
#include <vector>
#include <gsl/gsl_rng.h>			// random number generator gsl library        
#include <gsl/gsl_randist.h>		// gsl random distribution library
#include <math.h>					// standard math library			
using namespace	std;

// ________________________________________________________________________________________
// --------------------------------------------------------------------------Define classes

class TInd{					// Class defining a single individual
public:
	TInd();
	double lambda[2];											// mean fecundity
	double alpha[2];											// intraspecific competition coefficient
	double C_i[4][2];											// coefficient for the cubic calculating dispersal probability if infected 
	double C_s[4][2];										    // coefficient for the cubic calculating dispersal probability if suscpetible
	bool infectious_state;										// infection state
};

TInd::TInd(){				// Constructor for the class TInd
	lambda[0]=0;	
	lambda[1]=0;	
	alpha[0]=0;
	alpha[1]=0;
	for(int i=0;i<4;i++)
		for(int j=0;j<2;j++)
		{
			C_i[i][j]=0;
			C_s[i][j]=0;
		}
		infectious_state=0;
	}

class TPatch				// Class defining a single patch
{
public:
	TPatch();
	std::vector <TInd> females;									// vector of type TInd that stores all females in a patch
	std::vector <TInd> newfemales; 								// vector of type TInd that stores all new females in a patch during dispersal or reproduction
	std::vector <TInd> males;        							// vector of type TInd that stores all males in a patch
	std::vector <TInd> newmales;								// vector of type TInd that stores all new males in a patch during dispersal or reproduction
	double measured_dispersal;									// stores measured dipersal in a patch
	double measured_dispersal_infected;							// stores measured dispersal probability of infected individuals in a patch
	double measured_dispersal_susceptible;						// stores measured dispersal probability of susceptible individuals in a patch
	int No_of_infected;											// stores number of of infected individuals in a patch
	int No_of_susceptible;										// stores number of of susceptible individuals in a patch
	vector <double> mean_inf_allele;							// stores mean allelic values of cubic coefficients of infected individuals in a patch
	vector <double> mean_sus_allele;							// stores mean allelic values of cubic coefficients of susceptible individuals in a patch
};

TPatch::TPatch(){   		 // constructor for class TPatch
	females.clear();
	males.clear();
	newfemales.clear();
	newmales.clear();
	measured_dispersal=0;
	No_of_infected=0;
	No_of_susceptible=0;
	measured_dispersal_susceptible=0;
	measured_dispersal_infected=0;
	mean_inf_allele.clear();
	mean_sus_allele.clear();
}

// ________________________________________________________________________________________
// ----------------------------------------------------------------Declare global variables


int No;															// Initial number of individuals per patch
double LAMBDA;													// Mean fecundity
double ALPHA;													// Intraspecific competition coefficient
double DISPERSAL_PROB;											// Dispersal probability
int TMAX;														// Number of time steps
int REPLICATES;         										// Number of replicates
double DISP_MORT;												// dispersal mortality
double EXTINCTION_PROB;											// probabilty of local patch extinction
double VARIANCE; 												// variance of mutation
double MUT_RATE; 												// mutation rate for dispersal trait
double INIT_PREVALENCE; 										// Initial prevalence of disease
double TRANS_RATE;  											// transmission probability
double VIRULENCE; 												// fecundity virulence of parasite
double REL_MALE; 												// relative importance of male
double A; 														// search efficiency
bool CD; 														// condition dependency switch
bool GENETIC_VARIATION;											// switch for initial genetic variation
double HIST[21][2];												// histogram that stores the occurence of prevalences during the simulation
const int world_size = 10;										// number of patches in the metapopulation
TPatch world[world_size][world_size];							// Creating patches

const gsl_rng *gBaseRand;

// ________________________________________________________________________________________
// ------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);

	srand(randSeed);
	unsigned long r = rand();
	gsl_rng_set(gBaseRand, r);
}

// ________________________________________________________________________________________
// -------------------------------------------------------------------------Simplifications

// -------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

// ---------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

// -----------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}


const int RS = 100;                 							// random seed


//  random number generator
    // specify_rng(time(NULL));
    // specify_rng(RS);


// ________________________________________________________________________________________
// ------------------------------------------------------Functions for various calculations

double mean(double a, double b) 					 		// calculates the mean of two numbers
{
	return ((a+b)/double(2));
}

double mutate(double d, int t) 						 		// mutation procedure, returns modified allelic value
{
	double new_trait;
	double v=0;
	if(t <=TMAX*0.75)
		v=-((VARIANCE-0.1)*t/(TMAX*0.75))+VARIANCE;
	else v=0.1;
	if(ran()<MUT_RATE)
		new_trait= d+ gauss(v);
	else new_trait=d;
	return new_trait;

}

double densReg(double a) 					                // density regulation function
{
	return(1 /(1+(a)));
}

double infection_rate_calculator(double i,double s)	 		// calculates infection probability corresponding to a patch
{
	return 1-exp((-A*i*TRANS_RATE)/(1+(A*s)));
}

double effective_virulence(bool female, bool male)   		// calculates fecundity virulence
{
	double eff_vir;
	if(female==0 && male==0)
		eff_vir =double(1);
	if(female==1 && male==0)
		eff_vir = 1-VIRULENCE;
	if(female==0 && male==1)
		eff_vir = 1-(REL_MALE*VIRULENCE);
	if(female==1 && male==1)
		eff_vir = (1-(REL_MALE*VIRULENCE))*(1-VIRULENCE);
	return eff_vir;
}

void set_parameters()  										// read parameter input file  and set parameters
{
	string para[17];
	string line;
	ifstream myfile ("input.txt");							// open parameter input file
	int count=0;
	if (myfile.is_open())			
	{
		while ( getline (myfile,line))						// read the file
		{
			if(count%2==1)
			{
				para[count/2]=line;
			}
			count++;

		}
		myfile.close();
	}
	No = (int) std::atof(para[0].c_str());					// sets initial population density 
	LAMBDA=std::atof(para[1].c_str());						// sets mean fecundity
	ALPHA=std::atof(para[2].c_str());						// sets intraspecific competition coefficient
	TMAX= (int) std::atof(para[3].c_str());					// sets number of time steps fpr which the simulation runs
	REPLICATES=(int) std::atof(para[4].c_str());			// sets number of replicates
	DISPERSAL_PROB=std::atof(para[5].c_str());				// sets dispersal probability
	DISP_MORT=std::atof(para[6].c_str()); 					// sets dipersal mortality
	EXTINCTION_PROB=std::atof(para[7].c_str());				// sets probability of patch extinction
	VARIANCE=std::atof(para[8].c_str());					// sets standard deviation of gaussian for the mutation size
	MUT_RATE=std::atof(para[9].c_str());					// sets muation probability				
	INIT_PREVALENCE=std::atof(para[10].c_str());			// sets initial prevalence
	TRANS_RATE=std::atof(para[11].c_str());					// sets transmission rate
	VIRULENCE=std::atof(para[12].c_str());					// sets virulence
	REL_MALE=std::atof(para[13].c_str());					// sets relative contribution of male
	A=std::atof(para[14].c_str());							// sets search efficiency
	CD=(bool)std::atof(para[15].c_str());					// sets switch for condition dependency
	GENETIC_VARIATION=(bool)std::atof(para[16].c_str());	// sets switch for standing genetic variation
}


void init_world()			 		// initialise the metapopulation
{

	for(int x=0;x<world_size;x++)										 // clear all patches in the metapopulation
	{
		for(int y=0;y<world_size;y++)
		{
			world[x][y].females.clear();
			world[x][y].males.clear();
			for(int i=0; i<3;i++)
			{
				world[x][y].mean_inf_allele.push_back(0);
				world[x][y].mean_sus_allele.push_back(0);
			}
		}
	}
	for(int x=0;x<world_size;x++)
	{
		for(int y=0;y<world_size;y++)
		{

			for(int n=0; n<No; n++)										// creates new individual and assigns its attributes
			{
				TInd newind;							
				newind.lambda[0]=LAMBDA;
				newind.lambda[1]=LAMBDA;
				newind.alpha[0]=ALPHA;
				newind.alpha[1]=ALPHA;
				for(int i=0;i<3;i++)									
				{
					for(int j=0;j<2;j++)
					{
			  			if(GENETIC_VARIATION==1)						// if we are not initialising with standing genetic variation, all but C0 are set to zero
			  			{
			  				newind.C_s[i][j]= 0;
			  				newind.C_s[i][j]= 0;	
			  				if(CD==1)									// if dispersal depends on infection state
			  				{
			  					newind.C_i[i][j]=0; 	
			  					newind.C_i[i][j]=0; 	
			  				}
			  				else if(CD==0)								// if dispersal does not depend on infecton state
			  				{
			  					newind.C_i[i][j]= 33; 
			  					newind.C_i[i][j]= 33;	
			  				}
			  			}
			  			else if(GENETIC_VARIATION==0)					//if we are not initialsing with standing genetic variation set C0 to a fixed input					
			  			{
			  				newind.C_s[i][j]= DISPERSAL_PROB; 
			  				newind.C_s[i][j]= DISPERSAL_PROB;
			  				if(CD==1)
			  				{
			  					newind.C_i[i][j]=DISPERSAL_PROB; 
			  					newind.C_i[i][j]=DISPERSAL_PROB;
			  				}
			  				else if(CD==0)
			  				{
			  					newind.C_i[i][j]= 33; 
			  					newind.C_i[i][j]= 33;	
			  				}
			  			}
			  		}
			  	}
			  	for(int j=0;j<2;j++)									
			  	{
					if(GENETIC_VARIATION==1)							// if we are initialising with standing genetic, C0 is initialised from a uniform distribution between 0 and 1
					{
						newind.C_s[3][j]= ran(); 
						newind.C_s[3][j]= ran();	
						if(CD==1)
						{
							newind.C_i[3][j]=ran(); 	
							newind.C_i[3][j]=ran(); 	
						}
						else if(CD==0)
						{
							newind.C_i[3][j]= 33; 
							newind.C_i[3][j]= 33;	
						}
					}
					else if(GENETIC_VARIATION==0)
					{
						newind.C_s[3][j]= DISPERSAL_PROB; 
						newind.C_s[3][j]= DISPERSAL_PROB;
						if(CD==1)
						{
							newind.C_i[3][j]=DISPERSAL_PROB; 
							newind.C_i[3][j]=DISPERSAL_PROB;
						}
						else if(CD==0)
						{
							newind.C_i[3][j]= 33; 
							newind.C_i[3][j]= 33;	
						}
					}
				}
				if(ran()<INIT_PREVALENCE)
					newind.infectious_state=1;
				else newind.infectious_state=0;
				if(ran()<0.5)
					world[x][y].females.push_back(newind);   		// adding new females
				else world[x][y].males.push_back(newind);    		// adding new males
			}



		}
	}

	for(int h=0; h<=20; h++)
	{
		HIST[h][0]=double(h)/double(20);							// sets bins of the histogram of prevalences
		HIST[h][1]=0;												// initialises frequencies to zero
	}


}

void prevalence_measure()			// measures the number of susceptible and infected individuals in a patch
{

	for(int x=0;x<world_size;x++)
	{
		for(int y=0;y<world_size;y++)
		{
			int count_infected=0;
			int count_susceptible=0;
			for(int f=0;f<world[x][y].females.size();f++)
			{
				if(world[x][y].females.at(f).infectious_state ==1)
					count_infected++;
				else count_susceptible++;
			}
			for(int m=0;m<world[x][y].males.size();m++)
			{
				if(world[x][y].males.at(m).infectious_state ==1)
					count_infected++;
				else count_susceptible++;
			}

			world[x][y].No_of_susceptible=count_susceptible;
			world[x][y].No_of_infected=count_infected;
		}
	}
}

void hist()      				// makes a histogram of the prevalences in all patches throughout the simulation
{
	for(int x=0;x<world_size;x++)
	{
		for(int y=0; y<world_size;y++)
		{
			double prevalence=double(world[x][y].No_of_infected)/double(world[x][y].No_of_susceptible+world[x][y].No_of_infected);
			if(world[x][y].No_of_susceptible==0 && world[x][y].No_of_infected>0)
				HIST[19][1]=HIST[19][1]+1;
			else
			{
				for(int h=0;h<20;h++)
				{

					if(prevalence>=HIST[h][0] && prevalence<HIST[h+1][0])
						HIST[h][1]=HIST[h][1]+1;

				}
			}
		}
	}
	HIST[20][1]=HIST[19][1];
	
}

vector<int> decide_patch(int x, int y) 			// deciding coordinates of new patch-NN8 dispersal
{ 
	vector<int> coordinates;
	int newx=0;
	int newy=0;
	int decider=floor(ran()*double(8));			// draw 0,1,...,7 to choose one of the surrounding patches
	switch(decider){
		case 0: newx++; 						// patch to the right
		break;
		case 1: newy++;							// patch above
		break;
		case 2: newx++;							// patch diagonally above to the right
		newy++;
		break;
		case 3:newx--;							// patch to the left								
		break;
		case 4:newy--;							// patch below
		break;
		case 5:newx--;							// patch diagonally below to the left							
		newy--;
		break;
		case 6: newx++;							// patch diagonally below to the right
		newy--;
		break;
		case 7:newx--;							// patch diagonally above to the left
		newy++;
		break;
		default: cout<<"error in nn8"<<endl;
	}
	newx=newx+x;								// add the coordinate of the new patch
	newy=newy+y;
	if (newx<0)									// take care of the boundaries-periodic boundary condition
		newx=world_size - 1;
	if(newx==world_size)
		newx=0;
	if(newy<0)
		newy=world_size - 1;
	if(newy==world_size)
		newy=0;
	coordinates.push_back(newx);
	coordinates.push_back(newy);
	return coordinates;


}


void disperse()			                // dispersal procedure
{
	for(int x=0; x<world_size;x++)						// clears newmales and newfemales
	{
		for(int y=0; y<world_size;y++)
		{
			world[x][y].newfemales.clear();
			world[x][y].newmales.clear();
			world[x][y].mean_inf_allele.clear();
			world[x][y].mean_sus_allele.clear();
		}
	}
	for(int x=0;x<world_size;x++)						// dispersal loop
	{
		for(int y=0;y<world_size;y++)
		{
			int count_dispersers=0;
			int count_susceptible_disp=0;;
			int count_infected_disp=0;
			int count_ind=world[x][y].males.size()+world[x][y].females.size();
			double prevalence = double(world[x][y].No_of_infected)/double(world[x][y].No_of_infected + 
				world[x][y].No_of_susceptible);																													// calculates the prevalence in a patch
			vector <double> i_allele;
			vector<double> s_allele;
			for(int i=0; i<4;i++)
			{
				i_allele.push_back(0);
				s_allele.push_back(0);
			}
			for(int f=0;f<world[x][y].females.size();f++)
			{
				double dispersal_probability=0;
				for(int i=0; i<4;i++)
				{
					i_allele.at(i)=i_allele.at(i)+mean(world[x][y].females.at(f).C_i[i][0],
						world[x][y].females.at(f).C_i[i][1]);
					s_allele.at(i)=s_allele.at(i)+mean(world[x][y].females.at(f).C_s[i][0],
						world[x][y].females.at(f).C_s[i][1]);
				}
				if(CD==1)																																		
				{
					if(world[x][y].females.at(f).infectious_state==0)
						dispersal_probability=(mean(world[x][y].females.at(f).C_s[0][0],
							world[x][y].females.at(f).C_s[0][1])*prevalence*prevalence*prevalence)+(mean(world[x][y].females.at(f).C_s[1][0],					// calculate dispersal probability as a function of prevalence and infection state
							world[x][y].females.at(f).C_s[1][1])*prevalence*prevalence)+(mean(world[x][y].females.at(f).C_s[2][0],
							world[x][y].females.at(f).C_s[2][1])*prevalence)+mean(world[x][y].females.at(f).C_s[3][0],
							world[x][y].females.at(f).C_s[3][1]);
							else dispersal_probability=(mean(world[x][y].females.at(f).C_i[0][0],
								world[x][y].females.at(f).C_i[0][1])*prevalence*prevalence*prevalence)+(mean(world[x][y].females.at(f).C_i[1][0],
								world[x][y].females.at(f).C_i[1][1])*prevalence*prevalence)+(mean(world[x][y].females.at(f).C_i[2][0],
								world[x][y].females.at(f).C_i[2][1])*prevalence)+mean(world[x][y].females.at(f).C_i[3][0],
								world[x][y].females.at(f).C_i[3][1]);
							}
							else dispersal_probability=(mean(world[x][y].females.at(f).C_s[0][0],
								world[x][y].females.at(f).C_s[0][1])*prevalence*prevalence*prevalence)+(mean(world[x][y].females.at(f).C_s[1][0],
								world[x][y].females.at(f).C_s[1][1])*prevalence*prevalence)+(mean(world[x][y].females.at(f).C_s[2][0],
								world[x][y].females.at(f).C_s[2][1])*prevalence)+mean(world[x][y].females.at(f).C_s[3][0],
								world[x][y].females.at(f).C_s[3][1]);
								if(ran()< dispersal_probability)																									
								{
									if(world[x][y].females.at(f).infectious_state==0)																				
										count_susceptible_disp++;
									else count_infected_disp++;
								if(ran()>DISP_MORT)																												// if the individual does not die during dispersal put it in its target patch
								{
									std::vector<int> coor=decide_patch(x,y);
									world[coor.at(0)][coor.at(1)].newfemales.push_back(world[x][y].females.at(f));
								}
								world[x][y].females.erase(world[x][y].females.begin()+f);																		// erase the female which has either dispersed or is dead
								f--;																															// counter reduced to acount for female which has either dispersed or is dead
								count_dispersers++;
							}

						}
						for(int m=0;m<world[x][y].males.size();m++)
						{
							double dispersal_probability=0;
							for(int i=0; i<4;i++)
							{
								i_allele.at(i)=i_allele.at(i)+mean(world[x][y].males.at(m).C_i[i][0],world[x][y].males.at(m).C_i[i][1]);
								s_allele.at(i)=s_allele.at(i)+mean(world[x][y].males.at(m).C_s[i][0],world[x][y].males.at(m).C_s[i][1]);
							}
							if(CD==1)
							{
								if(world[x][y].males.at(m).infectious_state==0)
									dispersal_probability=(mean(world[x][y].males.at(m).C_s[0][0],
										world[x][y].males.at(m).C_s[0][1])*prevalence*prevalence*prevalence)+(mean(world[x][y].males.at(m).C_s[1][0],
										world[x][y].males.at(m).C_s[1][1])*prevalence*prevalence)+ (mean(world[x][y].males.at(m).C_s[2][0],
										world[x][y].males.at(m).C_s[2][1])*prevalence)+mean(world[x][y].males.at(m).C_s[3][0],
										world[x][y].males.at(m).C_s[3][1]);
										else dispersal_probability=(mean(world[x][y].males.at(m).C_i[0][0],
											world[x][y].males.at(m).C_i[0][1])*prevalence*prevalence*prevalence)+(mean(world[x][y].males.at(m).C_i[1][0],
											world[x][y].males.at(m).C_i[1][1])*prevalence*prevalence)+ (mean(world[x][y].males.at(m).C_i[2][0],
											world[x][y].males.at(m).C_i[2][1])*prevalence)+mean(world[x][y].males.at(m).C_i[3][0],
											world[x][y].males.at(m).C_i[3][1]);
										}
										else dispersal_probability=(mean(world[x][y].males.at(m).C_s[0][0],
											world[x][y].males.at(m).C_s[0][1])*prevalence*prevalence*prevalence)+(mean(world[x][y].males.at(m).C_s[1][0],
											world[x][y].males.at(m).C_s[1][1])*prevalence*prevalence)+ (mean(world[x][y].males.at(m).C_s[2][0],
											world[x][y].males.at(m).C_s[2][1])*prevalence)+mean(world[x][y].males.at(m).C_s[3][0],world[x][y].males.at(m).C_s[3][1]);
											if(ran()< dispersal_probability)
											{
												if(world[x][y].males.at(m).infectious_state==0)
													count_susceptible_disp++;
												else count_infected_disp++;
												if(ran()>DISP_MORT)
												{
													std::vector<int> coor=decide_patch(x,y);
													world[coor.at(0)][coor.at(1)].newmales.push_back(world[x][y].males.at(m));
												}
												world[x][y].males.erase(world[x][y].males.begin()+m);
												m--;
												count_dispersers++;
											}

										}

										world[x][y].measured_dispersal=double(count_dispersers)/double(count_ind);
										world[x][y].measured_dispersal_infected=double(count_infected_disp)/double(world[x][y].No_of_infected);
										world[x][y].measured_dispersal_susceptible=double(count_susceptible_disp)/double(world[x][y].No_of_susceptible);
										for(int i=0; i<4; i++)
										{
											world[x][y].mean_inf_allele.push_back(i_allele.at(i)/double(count_ind));
											world[x][y].mean_sus_allele.push_back(s_allele.at(i)/double(count_ind));
										}

									}

								}

								for(int x=0;x<world_size;x++)
								{
									for(int y=0;y<world_size;y++)
									{
										if(world[x][y].newfemales.size()>0)
										{
											for(int f=0;f<world[x][y].newfemales.size();f++)
											{
												world[x][y].females.push_back(world[x][y].newfemales.at(f));
											}
										}
										if(world[x][y].newmales.size()>0)
										{
											for(int m=0;m<world[x][y].newmales.size();m++)
											{
												world[x][y].males.push_back(world[x][y].newmales.at(m));
											}
										}
									}
								}

							}







void reproduce(int t) 			// reproduction loop
{ 
	for(int x=0; x<world_size;x++)
	{
		for(int y=0; y<world_size;y++)
		{
			world[x][y].newfemales.clear();
			world[x][y].newmales.clear();
			if(world[x][y].males.size()!=0 && world[x][y].females.size() != 0)
			{
				int Nf=world[x][y].females.size();
				int Nm=world[x][y].males.size();
				double alpha_net=0;
			  for(int f=0; f<Nf;f++)     																									 		// calculating net alpha fover all individuals
			  {
			  	alpha_net=alpha_net+mean(world[x][y].females.at(f).alpha[0],world[x][y].females.at(f).alpha[1]);
			  }
			  for(int m=0; m<Nm;m++)
			  {
			  	alpha_net=alpha_net+mean(world[x][y].males.at(m).alpha[0],world[x][y].males.at(m).alpha[1]);
			  }
			  for(int f=0;f<Nf;f++)
			  {
			  	int mate_position=floor(ran()*world[x][y].males.size());    																		// female chooses a mate 
			  	double mean_lambda= mean(world[x][y].females.at(f).lambda[0],world[x][y].females.at(f).lambda[1]);								
			      double mean_offspring = 2*mean_lambda * densReg(alpha_net)*effective_virulence(world[x][y].females.at(f).infectious_state,		// mean number of offspring after density regulation and fecundity virulence
			      	world[x][y].males.at(mate_position).infectious_state);											
			      int no_of_babies= poisson(mean_offspring);																						// number of offspring realised according to a poisson distribution
			      
			      for(int b=0;b<no_of_babies;b++)         																							// setting allelic values in the offspring, drawn at random, subjected to mutation
			      {
			      	TInd newind;
			      	newind.infectious_state=0;																										// all offspring are born susceptible
			      	if(ran()<0.5)
			      		newind.lambda[0]=world[x][y].females.at(f).lambda[0];																		// lambda and alpha are fixed, coefficients can mutate
			      	else newind.lambda[0]=world[x][y].females.at(f).lambda[1];
			      	if(ran()<0.5)
			      		newind.lambda[1]=world[x][y].males.at(mate_position).lambda[0];
			      	else newind.lambda[1]=world[x][y].males.at(mate_position).lambda[1];
			      	if(ran()<0.5)
			      		newind.alpha[0]=world[x][y].females.at(f).alpha[0];
			      	else newind.alpha[0]=world[x][y].females.at(f).alpha[1];
			      	if(ran()<0.5)
			      		newind.alpha[1]=world[x][y].males.at(mate_position).alpha[0];
			      	else newind.alpha[1]=world[x][y].males.at(mate_position).alpha[1];
			      	for(int i=0; i<4; i++)
			      	{
			      		if(CD==1)
			      		{
			      			if(ran()<0.5)
			      				newind.C_i[i][0]=mutate(world[x][y].females.at(f).C_i[i][0],t);
			      			else newind.C_i[i][0]=mutate(world[x][y].females.at(f).C_i[i][1],t);
			      			if(ran()<0.5)
			      				newind.C_i[i][1]=mutate(world[x][y].males.at(mate_position).C_i[i][0],t);
			      			else newind.C_i[i][1]=mutate(world[x][y].males.at(mate_position).C_i[i][1],t);
			      		}
			      		if(CD==0)
			      		{
			      			newind.C_i[i][1]=33;
			      			newind.C_i[i][0]=33;
			      		}

			      		if(ran()<0.5)
			      			newind.C_s[i][0]=mutate(world[x][y].females.at(f).C_s[i][0],t);
			      		else newind.C_s[i][0]=mutate(world[x][y].females.at(f).C_s[i][1],t);
			      		if(ran()<0.5)
			      			newind.C_s[i][1]=mutate(world[x][y].males.at(mate_position).C_s[i][0],t);
			      		else newind.C_s[i][1]=mutate(world[x][y].males.at(mate_position).C_s[i][1],t);
			      	}
			      	if(ran()<0.5)
			      		world[x][y].newfemales.push_back(newind);
			      	else world[x][y].newmales.push_back(newind);
			      }
			  }

			}
			else 
			{
				world[x][y].females.clear();
				world[x][y].males.clear();
			}


		}

	}

}

void transmission()								// transmission of parasite
{
	for(int x=0;x<world_size;x++)
	{
		for(int y=0;y<world_size;y++)
		{
			double infection_prob=  infection_rate_calculator(double(world[x][y].No_of_infected),			// calculation of infection probability accordinf to density of S and I in previous generation
				double(world[x][y].newmales.size()+world[x][y].newfemales.size()));
			for(int f=0;f<world[x][y].newfemales.size();f++)												// assign  infection states to all the new individuals
			{

				if(ran()<infection_prob)
					world[x][y].newfemales.at(f).infectious_state=1;

			}
			for(int m=0;m<world[x][y].newmales.size();m++)
			{
				if(ran()<infection_prob)
					world[x][y].newmales.at(m).infectious_state=1;

			}
		}
	}
}

void death()						// all individuals from previous generation die and are replaced by their offspring
{
	for(int x=0;x<world_size;x++)
	{
		for(int y=0;y<world_size;y++)
		{
			world[x][y].females.clear();
			world[x][y].males.clear();
			world[x][y].females=world[x][y].newfemales;
			world[x][y].males=world[x][y].newmales;
			world[x][y].newfemales.clear();
			world[x][y].newmales.clear();
		}
	}
}

void patch_extinction()			// local patch exinction with a given probability
{
	for(int x=0;x<world_size;x++)
	{
		for(int y=0; y<world_size;y++)
		{
			if(ran()<EXTINCTION_PROB)
			{
				world[x][y].females.clear();
				world[x][y].males.clear();
			}
		}
	}
}


int main()
{

	ofstream op;
	ofstream op1;
	ofstream op2;
	op.open("output.txt");																													// open output files
	op1.open("output_1.txt");
	op2.open("output_2.txt");
	specify_rng(RS);																														// set seed for random number generator
	set_parameters();
	op <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y"<<" "<<"N"<<" "<<"disp_rate"<<" "																// output the population size, number of S, I, measured dispersal, measured dispersal of S and I per patch
	<<"susceptible"<<" " <<"infected"<<" "<<"s_disp"<<" "<<"i_disp"<<" "<<endl;
	op1<<"rep"<<" "<<"s_0 i_0 s_1 i_1 s_2 i_2 s_3 i_3"<<endl;																				// output all coefficients for all the individuals at teh end of the simulation
	op2<<"rep"<<" "<<"prev"<<" "<<"freq"<<endl;																							    // output histogram
	op2<<endl;											
	for(int r=0; r<REPLICATES; r++)																											// replicates
	{
		init_world();
		for(int t=0; t<TMAX;t++)																											// time loop
		{
			prevalence_measure();	
			hist();
			if(t==TMAX-1)
			{
				for(int x=0; x<world_size;x++)
				{
					for(int y=0; y<world_size;y++)
					{
						op <<r<<" "<<t<<" "<<x<<" "<<y<<" "<<world[x][y].males.size()+world[x][y].females.size()<<
						" "<<world[x][y].measured_dispersal<< " "<<world[x][y].No_of_susceptible<<" "
						<<world[x][y].No_of_infected<<" "<<world[x][y].measured_dispersal_susceptible<<" "
						<<world[x][y].measured_dispersal_infected<<" ";
						op<<endl;
						for(int f=0;f<world[x][y].females.size();f++)
						{
							op1<<r<<" ";
							for(int i=0;i<4;i++)
							{
								op1<<mean(world[x][y].females.at(f).C_s[i][0],world[x][y].females.at(f).C_s[i][1])<<" ";
								op1<<mean(world[x][y].females.at(f).C_i[i][0],world[x][y].females.at(f).C_i[i][1])<<" ";
							}
							op1<<endl;
						}
						for(int m=0;m<world[x][y].males.size();m++)
						{
							op1<<r<<" ";
							for(int i=0;i<4;i++)
							{
								op1<<mean(world[x][y].males.at(m).C_s[i][0],world[x][y].males.at(m).C_s[i][1])<<" ";
								op1<<mean(world[x][y].males.at(m).C_i[i][0],world[x][y].males.at(m).C_i[i][1])<<" ";
							}
							op1<<endl;
						}

					}
				}
			}
			disperse();																													 // life cycle
			prevalence_measure();
			reproduce(t);
			transmission();
			death();
			patch_extinction();

		}


		for(int h=0;h<21;h++)
		{
			op2<<r<<" ";
			op2<<HIST[h][0]<<" "<<HIST[h][1];
			op2<<endl;
		}
	}
	

	op.close();																															// close output files
	op1.close();
	op2.close();
	return 0;

}




