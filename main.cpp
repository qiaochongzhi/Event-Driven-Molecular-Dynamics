/* This program,we calculate the thermodynamic properties of
a confined fluid in a random porous media or slit pore by 
Molecular dynamic simulation. */

// v 0.1.0

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <math.h>
#include "Randnumber.h"

#include <sstream>
/************************************ Period parameters ***********************************/
#define FALSE 0
#define TRUE 1
#define MAX_NUMBER_OF_PARTICLES 10000 /*The maximum number of the particles of each species*/
#define PI 3.14159265
#define BLOCK 1000
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define Length sizeof(struct DATA)
#define INFINT 1.e10
#define BOUNDARY 300.0
#define BIGTIME 300.0

using namespace std;

string int_to_string(int i)
{
  char ch[10];
  sprintf(ch, "%d", i);
  string s(ch);

  return s;
}

typedef struct
{
  double x;
  double y;
  double z;
} VECTOR;

typedef struct
{
  double time;
  int i; 
  /* -1 is the left wall;
   * -2 is the right wall;
   * -3 is the left square well;
   * -4 is the right square well.
  */
} RESULT_WALL;

typedef struct DATA
{
  double x;
  struct DATA * next;
};

class particle
{
public:
  void setParticle(double, int);
  particle() {}
  ~particle() 
  {
    delete[]position;
    delete[]velocity; 
    delete[]collisionTime;
    delete[]partner;
  }
  string type;

  double diameter;
  int number;
  VECTOR *velocity;
  VECTOR *position;
  double *collisionTime;
  int *partner;
};

class diffusionBuffer
{
public:
  void initialDiffusionBuffer(int number)
  {
    this->count = 0;
    this->originalPosition = new VECTOR[number];
    this->truePosition = new VECTOR[number];
    this->distance = new double[number];

    for (int i = 0; i < number; i++)
    {
      originalPosition[i].x = 0;
      originalPosition[i].y = 0;
      originalPosition[i].z = 0;

      truePosition[i].x = 0;
      truePosition[i].y = 0;
      truePosition[i].z = 0;

      distance[i] = 0;
    }
  }
  diffusionBuffer() {}
  ~diffusionBuffer()
  {
    delete[]originalPosition;
    delete[]truePosition;
    delete[]distance;
  }

  int count;
  VECTOR *originalPosition;
  VECTOR *truePosition;
  double *distance;
};

void particle::setParticle(double diameter, int number)
{
  this->diameter = diameter;
  this->number = number;
  this->position = new VECTOR[number];
  this->velocity = new VECTOR[number];
  this->partner = new int[number];
  this->collisionTime = new double[number];
}

class MD
{
public:
  MD();
  ~MD();
  void MDrun();
  void printPoreDensityDistributionTemp(int);
private:
  int isOverlap(VECTOR, VECTOR, double, double);
  void boundaryConditions(VECTOR* dr)
  {
    if (dr->x > 0.5*L.x)  dr->x -= L.x;
    if (dr->x < -0.5*L.x) dr->x += L.x;

    if (dr->y > 0.5*L.y)  dr->y -= L.y;
    if (dr->y < -0.5*L.y) dr->y += L.y;

    if (!isSlitPore)
    {
      if (dr->z > 0.5*L.z)  dr->z -= L.z;
      if (dr->z < -0.5*L.z) dr->z += L.z;
    }
  }
  //void dispalceParticles(int);
  void initialParticles();
  void initialCollisionTimeList();
  void updateCollisionTimeList(int);
  void forwardTime(double);
  double bump(int);
  double calculateCollision(int, int);

  void poreDensityDistribution(double*);
  void densityDistribution(double*);
  void printPoreDensityDistribution();
  void printDensityDestribution();

  void evaluteDiffusion();
  void initialDiffusion();
  void accumulateDiffusion();
  void printDiffusion();
  void zeroDiffusion();

  RESULT_WALL calculateCollision(int);

  double numberOfCycles;
  double numberOfInitializationCycles;
  double particleForm;
  double beta;

  double deltaTime;
  //double limitDiffuseTime;
  int numberOfBuffer;
  int numberOfDiffusionStep; 
  int countAvarageOfDiffusion;

  // not initial
  double* averageDiffusionDistance;
  int limitAverageOfDiffusion;
  int diffusionStep;

  double* density;
  double* density1;
  double unitSize;
  int numberOfNode;

  double chemP;
  double errorChemp;

  int isSlitPore;
  double depthOfWell;
  double widthOfWell;

  string printForm;
  string slitPoreForm;

  VECTOR L;

  particle particles;
  diffusionBuffer* buffer;

  clock_t timeEnd;
  clock_t timeStart;

};

MD::MD()
{
  char tip[50];
  double diameter;
  string s;
  string ss;
  int number;
  
  timeStart = clock();

  cout << "start" << endl;

  ifstream input("input.dat", ios::in);

  getline(input, s);
  input >> L.x >> L.y >> L.z;
  getline(input, s);

  getline(input, s);
  input >> numberOfNode;  
  
  input >> tip;
  input >> isSlitPore;

  density = new double[numberOfNode];
  density1 = new double[numberOfNode];
  if (isSlitPore)
  {
    for (int i = 0; i < numberOfNode; i++)
    {
      density[i] = 0;
    }
  } 
  else
  {
    for (int i = 0; i < numberOfNode; i++)
    {
      density1[i] = 0;
    }

  }

  if (isSlitPore)
  {
    unitSize = L.z / numberOfNode;
  } 
  else
  {
    unitSize = L.z / (numberOfNode * 2);
  }
 

  
  //cout << isSlitPore << endl;

  input >> tip;
  input >> tip >> tip;
  input >> diameter >> number;
  getline(input, s);

  particles.setParticle(diameter, number);

  //chemP = new double[number];
  //errorChemp = new double[number];

  getline(input, s);
  input >> particleForm;
  getline(input, s);

  getline(input, s);
  input >> slitPoreForm;
  getline(input, s);

  getline(input, s);
  input >> depthOfWell >> widthOfWell;

  input >> tip;
  input >> numberOfCycles;
  numberOfCycles *= BLOCK;

  input >> tip;
  input >> numberOfInitializationCycles;
  numberOfInitializationCycles *= BLOCK;

  input >> tip;
  input >> beta;

  getline(input, s);
  getline(input, s);
  input >> printForm;

  getline(input, s);
  getline(input, s);
  input >> deltaTime;

  getline(input, s);
  getline(input, s);
  input >> numberOfDiffusionStep;
  numberOfDiffusionStep *= BLOCK;

  averageDiffusionDistance = new double[numberOfDiffusionStep];
  for (int i = 0; i < numberOfDiffusionStep; i++)
    averageDiffusionDistance[i] = 0;

  getline(input, s);
  getline(input, s);
  input >> numberOfBuffer;

  if (numberOfBuffer > 0)
  {
    buffer = new diffusionBuffer[numberOfBuffer];
    for (int i = 0; i < numberOfBuffer; i++)
      buffer[i].initialDiffusionBuffer(numberOfDiffusionStep);
  }

  getline(input, s);
  getline(input, s);
  input >> limitAverageOfDiffusion;

  getline(input, s);
  getline(input, s);
  input >> diffusionStep;

  initialDiffusion();
  
  /*{
    cout << "limit average of diffusion: " << limitAverageOfDiffusion << endl;
    exit(1);
  }*/
  /*{
    cout << depthOfWell << widthOfWell << endl;
    cout << beta << endl;
  } */

  input.close();
}

MD::~MD()
{
  delete[]buffer;
  delete[]averageDiffusionDistance;

  //delete[]chemP;
  //delete[]errorChemp;

  if (isSlitPore)
  {
    delete[]density;
  } 
  else
  {
    delete[]density1;
  }
}

int MD::isOverlap(VECTOR p1, VECTOR p2, double sigma1, double sigma2)
{
  int overlap;
  VECTOR dr;

  dr.x = p1.x - p2.x;
  dr.y = p1.y - p2.y;
  dr.z = p1.z - p2.z;

  boundaryConditions(&dr);

  double d = SQR(dr.x) + SQR(dr.y) + SQR(dr.z) - SQR((sigma1 + sigma2) / 2);

  overlap = d < 0 ? TRUE : FALSE;

  return overlap;
}

void MD::initialParticles()
{
  double packingFraction = 0;
  double sigma;
  double totalEnergy = 0.0;
  double energyFactor;
  VECTOR trial;
  int overlap;
  
  packingFraction += PI*CUBE(particles.diameter)*particles.number / 6;

  packingFraction /= L.x * L.y * L.z;
  if (packingFraction > 1.0)
  {
    cout << "Error in the initial number of particles" << endl;
    exit(-1);
  }

  for (int j = 0;j < particles.number;j++)
  {
    trial.x = (RandomNumber() - 0.5)*L.x;
    trial.y = (RandomNumber() - 0.5)*L.y;
    if (isSlitPore)
    {
      trial.z = (RandomNumber() - 0.5)*(L.z - particles.diameter);
      //trial.z = -RandomNumber()*(L.z - particles.diameter)/2;
    } 
    else
    {
      trial.z = (RandomNumber() - 0.5)*L.z;
    }

    sigma = particles.diameter;

    overlap = FALSE;

    for (int l = 0;l < j&&overlap == FALSE;l++)
    {
      overlap = isOverlap(trial, particles.position[l], sigma, sigma);
    }

    if (overlap == FALSE)
    {
      particles.position[j].x = trial.x;
      particles.position[j].y = trial.y;
      particles.position[j].z = trial.z;

      particles.velocity[j].x = RandomNumber() - 0.5;
      particles.velocity[j].y = RandomNumber() - 0.5;
      particles.velocity[j].z = RandomNumber() - 0.5;

      totalEnergy += SQR(particles.velocity[j].x);
      totalEnergy += SQR(particles.velocity[j].y);
      totalEnergy += SQR(particles.velocity[j].z);
    }
    else
    {
      j--;
    }
  }
  energyFactor = sqrt(3 * particles.number / (beta * totalEnergy));
  for (int i = 0; i < particles.number; i++) 
  {
    particles.velocity[i].x *= energyFactor;
    particles.velocity[i].y *= energyFactor;
    particles.velocity[i].z *= energyFactor;

    particles.collisionTime[i] = BIGTIME;
    particles.partner[i] = i;
  }

  /*{
    double Ttest = 0;
    for (int i = 0; i < particles.number; i++)
    {
      Ttest += SQR(particles.velocity[i].x);
      Ttest += SQR(particles.velocity[i].y);
      Ttest += SQR(particles.velocity[i].z);
    }
    Ttest /= 3 * particles.number;
    cout << "Ttest = " << Ttest << endl;
  } */

}

void MD::initialCollisionTimeList()
{
//  VECTOR dr;
//  VECTOR dv;
//  double BIJ;
//  double RijSqr;
//  double VijSqr;
//  double disCR;
  double tij;
  RESULT_WALL tw;

  for (int i = 0; i < particles.number-1; i++) 
  {
    tw = calculateCollision(i);
    if (tw.time < particles.collisionTime[i])
    {
      particles.collisionTime[i] = tw.time;
      particles.partner[i] = tw.i;
    }
    for (int j = i; j < particles.number; j++) 
    {
      tij = calculateCollision(i, j);
      if (tij < particles.collisionTime[i])
      {
        particles.collisionTime[i] = tij;
        particles.partner[i] = j;
      }
    }
  }
}

double MD::calculateCollision(int i, int j)
{
  VECTOR dr;
  VECTOR dv;
  double BIJ;
  double RijSqr;
  double VijSqr;
  double disCR;
  double tij = BIGTIME;

  dr.x = particles.position[i].x - particles.position[j].x;
  dr.y = particles.position[i].y - particles.position[j].y;
  dr.z = particles.position[i].z - particles.position[j].z;
  boundaryConditions(&dr);

  dv.x = particles.velocity[i].x - particles.velocity[j].x;
  dv.y = particles.velocity[i].y - particles.velocity[j].y;
  dv.z = particles.velocity[i].z - particles.velocity[j].z;

  BIJ = dr.x * dv.x + dr.y * dv.y + dr.z * dv.z;

  if (BIJ < 0) 
  {
    RijSqr = SQR(dr.x) + SQR(dr.y) + SQR(dr.z);
    VijSqr = SQR(dv.x) + SQR(dv.y) + SQR(dv.z);
    disCR = SQR(BIJ) - VijSqr * (RijSqr - SQR(particles.diameter));
    if (disCR > 0)
    {
      tij = (- BIJ - sqrt(disCR)) / VijSqr;
    }
  }
  return tij;
}

RESULT_WALL MD::calculateCollision(int i) 
/* TODO: This stupid part should be modified. */
{
  RESULT_WALL tw;

  double rz;
  double vz;

  //double t1, t2;

  vz = particles.velocity[i].z;
  rz = particles.position[i].z;

  if (vz < 0)
  {
    if (rz > L.z/2 - widthOfWell && depthOfWell != 0)
    {
      tw.time = abs(rz - (L.z / 2 - widthOfWell)) / abs(vz);
      tw.i = -4; // exit from right square well;
    } 
    else if (rz > -(L.z / 2 - widthOfWell) && depthOfWell != 0)
    {
      tw.time = abs(rz + (L.z / 2 - widthOfWell)) / abs(vz);
      tw.i = -3; // enter the left square well;
    }
    else
    {
      tw.time = abs(rz + (L.z - particles.diameter) / 2) / abs(vz);
      tw.i = -1; // collision with left wall;
    }
  }
  else
  {
    if (rz < -(L.z / 2 - widthOfWell) && depthOfWell != 0)
    {
      tw.time = abs(rz + (L.z / 2 - widthOfWell)) / abs(vz);
      tw.i = -3; // exit from left square well;
    }
    else if (rz < (L.z / 2 - widthOfWell) && depthOfWell != 0)
    {
      tw.time = abs(rz - (L.z / 2 - widthOfWell)) / abs(vz);
      tw.i = -4; // enter the right square well;
    }
    else
    {
      tw.time = abs(rz - (L.z-particles.diameter) / 2) / abs(vz);
      tw.i = -2; // collision with right wall;
    }
  }

  return tw;
}

void MD::forwardTime(double dt)
{
  for(int i = 0; i < particles.number; i++)
  {
    particles.collisionTime[i] -= dt;
    particles.position[i].x += particles.velocity[i].x * dt; 
    particles.position[i].y += particles.velocity[i].y * dt; 
    particles.position[i].z += particles.velocity[i].z * dt; 

    boundaryConditions(&particles.position[i]);
  }
}

void MD::updateCollisionTimeList(int i)
{
  int ii = i;
  int jj = particles.partner[i];
  int rule;
  int init;
  double tij;

  RESULT_WALL tw;

  for (int i = 0; i < particles.number; i++) 
  {
    int k = particles.partner[i];
    rule = (i == ii)||(i == jj);
    rule = rule||(k == ii);
    rule = rule||(k == jj);
    if (rule) 
    {
      particles.collisionTime[i] = BIGTIME;
      if (isSlitPore)
      {
        tw = calculateCollision(i);
        if (tw.time < particles.collisionTime[i])
        {
          particles.collisionTime[i] = tw.time;
          particles.partner[i] = tw.i;
        }
      }
      if (i == ii || i == jj)
      {
        init = 0;
      }
      else 
      {
        init = i + 1; 
      }
      for (int j = init; j <particles.number; j++)
      {
        if (j != i)
        {
          tij = calculateCollision(i,j);
          
          if (tij < particles.collisionTime[i])
          {
            particles.collisionTime[i] = tij;
            particles.partner[i] = j;
          }
          if (tij < particles.collisionTime[j])
          {
            particles.collisionTime[j] = tij;
            particles.partner[j] = i;
          }
        }
      }
    }
  }

}

double MD::bump(int i)
{
  VECTOR dr;
  VECTOR dv;
  VECTOR dvv;
  int j = particles.partner[i];
  double factor;
  double sigma = particles.diameter;
  double b;
  double w;

  if (j >= 0)
  {
    dr.x = particles.position[i].x - particles.position[j].x;
    dr.y = particles.position[i].y - particles.position[j].y;
    dr.z = particles.position[i].z - particles.position[j].z;
    boundaryConditions(&dr);

    dv.x = particles.velocity[i].x - particles.velocity[j].x;
    dv.y = particles.velocity[i].y - particles.velocity[j].y;
    dv.z = particles.velocity[i].z - particles.velocity[j].z;

    factor = (dr.x * dv.x + dr.y * dv.y + dr.z * dv.z) / SQR(sigma);

    dvv.x = -factor * dr.x;
    dvv.y = -factor * dr.y;
    dvv.z = -factor * dr.z;

    particles.velocity[i].x += dvv.x;
    particles.velocity[j].x -= dvv.x;

    particles.velocity[i].y += dvv.y;
    particles.velocity[j].y -= dvv.y;

    particles.velocity[i].z += dvv.z;
    particles.velocity[j].z -= dvv.z;

    w = dvv.x * dr.x + dvv.y * dr.y + dvv.z * dr.z;
  } 
  else
  {
    w = 0;

    if (j == -1)
    {
      particles.velocity[i].z = -particles.velocity[i].z;
    }
    else if (j == -2)
    {
      particles.velocity[i].z = -particles.velocity[i].z;
    }
    else if (j == -3)
    {
      if (particles.velocity[i].z > 0)
      {
        // exit from the left square well;
        b = SQR(particles.velocity[i].z) + (2 * depthOfWell);
        if (b < 0)
        {
          particles.velocity[i].z = -particles.velocity[i].z;
        } 
        else
        {
          particles.velocity[i].z = sqrt(b);
        }
      } 
      else
      {
        // enter the left square well;
        b = SQR(particles.velocity[i].z) - (2 * depthOfWell);
        if (b < 0)
        {
          particles.velocity[i].z = -particles.velocity[i].z;
        } 
        else
        {
          particles.velocity[i].z = -sqrt(b);
        }
      }
    }
    else if (j == -4)
    {
      if (particles.velocity[i].z > 0)
      {
        // enter the right square well;
        b = SQR(particles.velocity[i].z) - (2 * depthOfWell);
        if (b < 0)
        {
          particles.velocity[i].z = -particles.velocity[i].z;
        } 
        else
        {
          //cout << "vz = " << particles.velocity[i].z << endl;
          //cout << "new vz = " << sqrt(b) << endl;
          particles.velocity[i].z = sqrt(b);
          
        }
      }
      else
      {
        // exit the right square well;
        b = SQR(particles.velocity[i].z) + (2 * depthOfWell);
        if (b < 0)
        {
          particles.velocity[i].z = -particles.velocity[i].z;
        } 
        else
        {
          //cout << "vz = " << particles.velocity[i].z << endl;
          //cout << "new vz = " << sqrt(b) << endl;
          particles.velocity[i].z = -sqrt(b);
        }
      }
    }
    else
    {
      cout << "error: [line 556] unknown collection type." << endl;
      exit(-1);
    }
  }

  return w;
}

void MD::MDrun()
{
  double totalTime = 0;
  double tij = BIGTIME;
  double PVNKT1;
  int nextCollisionParticle;
  double acw = 0;
  double tnow = 0;

  double nCount = 0;
  
  // Chemical potential //
  double chemPSum;
  double averageChemP;
  double enn;

  DATA *p;
  DATA *pHead;
  DATA *pTemp;

  p = new DATA;
  pHead = new DATA;
  pTemp = new DATA;

  chemP = 0;
  averageChemP = 0;
  chemPSum = 0;

  p = (struct DATA*)malloc(Length);
  pHead = p;

  VECTOR testParticle;
  // chemocal potential end //

  initialParticles();
  initialCollisionTimeList();

  cout << "begin" << endl;
  
  for (int i = 0; i < numberOfCycles; i++) 
  {
    tij = BIGTIME;
    if (i%1000 == 0)
    {
      cout << i << endl;

    }
    //if (i%1000 == 0)
    //{
    //  printPoreDensityDistributionTemp(i/1000);
    //}

    tnow = 0;
    
    while (1)
    {
      tij = BIGTIME;
      for (int j = 0; j < particles.number; j++)
      {
        if (particles.collisionTime[j] < tij)
        {
          tij = particles.collisionTime[j];
          nextCollisionParticle = j;
        }
      }

      if (tnow + tij > deltaTime)
      {
        forwardTime(deltaTime - tnow);
        //cout << "time now: " << tnow << endl;
        break;
      }

      //if (nextCollisionParticle < 0)
      //  cout << nextCollisionParticle << endl;

      totalTime += tij;
      forwardTime(tij);

      tnow += tij;
      //cout << "delta time: " << dt << endl;

      acw += bump(nextCollisionParticle);
      //cout << "new" << endl;
      //cout << particles.position[nextCollisionParticle].z << endl;
      //cout << nextCollisionParticle << endl;
      //cout << acw << endl;
      updateCollisionTimeList(nextCollisionParticle);
    
    }

    if (i > numberOfInitializationCycles)
    {
      nCount++;
      if (isSlitPore)
      {
        poreDensityDistribution(density);
      }
      else
      {
        densityDistribution(density1);
      }
      //cout << i / diffusionStep << endl;
      if (i%diffusionStep == 0)
      {
        evaluteDiffusion();
      }

      // chemcial potential //
      {
        testParticle.x = (RandomNumber() - 0.5)*L.x;
        testParticle.y = (RandomNumber() - 0.5)*L.y;
        testParticle.z = (RandomNumber() - 0.5)*L.z;

        enn = 1;
        for (int i = 0; i < particles.number; i++)
        {
          if (isOverlap(testParticle, particles.position[i], particles.diameter, particles.diameter))
          {
            enn = 0;
            break;
          }
        }


        chemPSum += enn;
        p->x = enn;
        pTemp = p;
        p = (struct DATA*)malloc(Length);
        pTemp->next = p;
      }
      // chemical potential end //

    }
    
  }

  // chemical potential //
  pTemp->next = NULL;
  free(p);
  averageChemP = chemPSum / nCount;
  chemP = -log(averageChemP);

  p = pHead;
  // counter = 0;
  // countBlock = 0;
  // sumBlock = 0;

  cout << "Chemical potential: " << chemP << endl;

  // end //


  PVNKT1 = beta * acw / ( (double)particles.number * 3.0 * totalTime);
  
  cout << "PV/NKT - 1 is : " << PVNKT1 << endl;
  if (isSlitPore)
  {
    for (int i = 0; i < numberOfNode; i++)
    {
      density[i] = density[i] / (nCount*unitSize*L.x*L.y);
    }
  }
  else
  {
    double dis = 0;
    for (int i = 0; i < numberOfNode; i++)
    {
      dis = (4 / 3)*PI*pow(unitSize, 3) *(pow(i + 1, 3) - pow(i, 3));
      density1[i] = density1[i] / (nCount*dis);
    }
  }
  if (isSlitPore)
  {
    printPoreDensityDistribution();
  }
  else
  {
    printDensityDestribution();
  }

  timeEnd = clock();
}

void MD::printPoreDensityDistributionTemp(int kk)
{
  int num = particles.number;
  int l;
  for (int i = 0; i < numberOfNode; i++)
  {
    density[i] = 0;
  }

  string s, sout1, sout2;
  sout1 = "Temp/density";
  sout2 = ".dat";

  s = int_to_string(kk + 1);
  s = sout1 + s + sout2;


  for (int i = 0; i < num; i++)
  {
    l = floor((particles.position[i].z + L.z / 2) / unitSize);
    //cout << poreSize << "  " << unitSize << "    " << particles[k].position[i].z << "   " << l << endl;
    density[l] += 1;
  }
  //for (int i = 0; i < numberOfNode; i++)
  //{
  //  cout << density[i] << endl;
  //}
  
  ofstream op2(s, ios::out);

  op2 << "Z" << '\t' << "fluid" << endl;

  for (int i = 0; i < numberOfNode; i++)
  {
    //op2 << (i + 0.5)*unitSize;
    op2 << (i)*unitSize << '\t' << density[i] / (unitSize*L.x*L.y) << endl;
  }
  op2.close();

}

void MD::poreDensityDistribution(double* d)
{
  int num = particles.number;
  int l;

  for (int i = 0; i < num; i++)
  {
    l = floor((particles.position[i].z + L.z / 2) / unitSize);
    //cout << poreSize << "  " << unitSize << "    " << particles[k].position[i].z << "   " << l << endl;
    d[l] += 1;
  }
}

void MD::printPoreDensityDistribution()
{
  ofstream op1("density.dat", ios::out);

  op1 << "Z" << '\t' << "fluid" << endl;

  for (int i = 0; i < numberOfNode; i++)
  {
    //op1 << (i + 0.5)*unitSize;
    op1 << (i)*unitSize << '\t' << density[i] << endl;
  }
  op1.close();
}

void MD::densityDistribution(double* d)
{
  int l;
  double dis;
  VECTOR t, trial;

  t.x = particles.position[0].x;
  t.y = particles.position[0].y;
  t.z = particles.position[0].z;

  for (int j = 0; j < particles.number; j++)
  {
    if (j == 0)
      continue;
    trial.x = t.x - particles.position[j].x;
    trial.y = t.y - particles.position[j].y;
    trial.z = t.z - particles.position[j].z;

    boundaryConditions(&trial);

    dis = sqrt(SQR(trial.x) + SQR(trial.y) + SQR(trial.z));

    if (dis < L.z / 2)
    {
      l = floor(dis / unitSize);
      if (l < numberOfNode)
      {
        density1[l] += 1;
        //cout << l << '\t' << density1[l] << endl;
      }
    }
  }
}

void MD::printDensityDestribution()
{
  string s, sout1, sout2;
  sout1 = "density";
  sout2 = ".dat";

  {
    //s = int_to_string(i + 1);
    //s = sout1 + s + sout2;

    ofstream op1("density1.dat", ios::out);

    op1 << "Z" << '\t' << "fluid" << endl;
    for (int j = 0; j < numberOfNode; j++)
    {
      op1 << (j + 0.5)*unitSize << '\t' << density1[j] << endl;
    }
    op1.close();
  }
}

void MD::initialDiffusion()
{
  for (int nb = 0; nb < numberOfBuffer; nb++)
  {
    buffer[nb].count = -nb * numberOfDiffusionStep / numberOfBuffer;
    zeroDiffusion();
  }
}

void MD::zeroDiffusion()
{
  countAvarageOfDiffusion = 0;
  for (int i = 0; i < numberOfDiffusionStep; i++)
    averageDiffusionDistance[i] = 0;
}

void MD::evaluteDiffusion()
{
  VECTOR dr;
  int ni;

  for (int nb = 0; nb < numberOfBuffer; nb++)
  {
    if (buffer[nb].count == 0)
    {
      for (int n = 0; n < particles.number; n++)
      {
        buffer[nb].originalPosition[n] = particles.position[n];
        buffer[nb].truePosition[n] = particles.position[n];
      }
    }
    if (buffer[nb].count >= 0)
    {
      ni = buffer[nb].count;
      buffer[nb].distance[ni] = 0.0;

      for (int n = 0; n < particles.number; n++)
      {
        dr.x = buffer[nb].truePosition[n].x - particles.position[n].x;
        dr.y = buffer[nb].truePosition[n].y - particles.position[n].y;
        dr.z = buffer[nb].truePosition[n].z - particles.position[n].z;

        dr.x = L.x * round(dr.x / L.x);
        dr.y = L.y * round(dr.y / L.y);
        dr.z = L.z * round(dr.z / L.z);

        buffer[nb].truePosition[n].x = particles.position[n].x + dr.x;
        buffer[nb].truePosition[n].y = particles.position[n].y + dr.y;
        buffer[nb].truePosition[n].z = particles.position[n].z + dr.z;

        dr.x = buffer[nb].truePosition[n].x - buffer[nb].originalPosition[n].x;
        dr.y = buffer[nb].truePosition[n].y - buffer[nb].originalPosition[n].y;
        dr.z = buffer[nb].truePosition[n].z - buffer[nb].originalPosition[n].z;

        buffer[nb].distance[ni] += (SQR(dr.x) + SQR(dr.y) + SQR(dr.z));
      }
    }
    buffer[nb].count++;
  }
  accumulateDiffusion();
}

void MD::accumulateDiffusion()
{
  double factor;
  
  for (int nb = 0; nb < numberOfBuffer; nb++)
  {
    if (buffer[nb].count == numberOfDiffusionStep)
    {
      for (int j = 0; j < numberOfDiffusionStep; j++)
      {
        averageDiffusionDistance[j] += buffer[nb].distance[j];
      }
      buffer[nb].count = 0;
      countAvarageOfDiffusion++;
      cout << countAvarageOfDiffusion << endl;
      if (countAvarageOfDiffusion == limitAverageOfDiffusion)
      {
        factor = 1.0 / (6 * numberOfDiffusionStep * particles.number * diffusionStep * deltaTime * limitAverageOfDiffusion);
        for (int j = 0; j < numberOfDiffusionStep; j++)
        {
          averageDiffusionDistance[j] *= factor ;
        }
        printDiffusion();
        zeroDiffusion();
      }
    }
  }
}

void MD::printDiffusion()
{
  //string s;

  {
    ofstream op1("duffision.dat", ios::out);

    op1 << "t" << '\t' << "distance" << endl;
    for (int j = 0; j < numberOfDiffusionStep; j++)
    {
      op1 << (j * deltaTime * diffusionStep) << '\t' << averageDiffusionDistance[j] << endl;
    }
    op1.close();
  }
}
  

int main()
{
  MD test;
  test.MDrun();

  cout << "done" << endl;
  return 0;
}
