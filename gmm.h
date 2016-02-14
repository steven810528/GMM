#ifndef GMM_H
#define GMM_H

#include <math.h>
#include <cmath>
#include <ctime>

class GMM
{
    public:
        vector<string> log;
    
        int cluster;
        const static int dim = 12;
        float convergence;
        //static int counter;
        vector<float> prePara;
    
        vector<vector<float> > nowBeta;
        vector< vector<float> > mean;
        vector<float> var; //segma^2
        vector<float> part;
    
        vector<vector<float> > data;
    
        GMM(int cluster, vector<vector<float> > &data);
        void init(int cluster);
        void calculate();
        void printPar();
        void save();
        void load();
    
    //private:
        void calBeta();
        void updatePart();
        void updateMean();
        void updateVar();
        void Kmean(float std);
        float diff();
        void normalize();
        float beta(int f,int c);
        float density(int f,int c);
        float distance(int f,int c);
        float distance(int f,vector<float> m );
};


#endif

