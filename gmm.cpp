#include "gmm.h"
#include "../common/print.h"
#include "../common/toolbox.h"
/*

static int cluster;
static float convergence;
vector<float> part;
vector< vector<float> > mean;
vector<float> var;
vector<vector<float> > *data;*/

time_t now; // 變數宣告
void printLog(vector<string> v);

void multipleVector(vector<float>& v, float m,int dim);
void addVector(vector<float>& v1 ,vector<float>v2 ,int dim );
void divVector(vector<float>& v1 ,float total ,int dim );
//float distance(int x,int c);
//float gauss_distribution(float mean, float std);

GMM::GMM(int cluster, vector<vector<float> > &v)
{
    //init
    //=1000;
    this->convergence=0.01;
    this->data=v;
    this->normalize();
    init(cluster);
    //printBigVector(this->data);
}
void GMM::init(int cluster)
{
    //宣告空間
    this->cluster=cluster;
    this->part.resize(cluster);
    this->mean.resize(cluster);
    for (int i=0;i<cluster;i++)
    {
        this->mean[i].resize(dim);
    }
    this->var.resize(cluster);
    
    nowBeta.resize(data.size());
    for(int i=0;i<nowBeta.size();i++)
    {
        nowBeta[i].resize(cluster);
    }
    //設定初始值
    for(int i =0;i<cluster;i++)
    {
        this->part[i]=(float)1/cluster;
        //this->var[i]=5000;
    }
    this->Kmean(2);
    /*
    //設定mean
    vector<float> tmp;
    tmp.resize(dim);
    for(int i=0;i<data.size();i++)
    {
        for(int j=0;j<dim;j++)
        {
            tmp[j]+=data[i][j];
        }
    }
    //cout<<tmp[0]<<endl;
    for(int c=0;c<cluster;c++)
    {
        for(int i=0;i<dim;i++)
        {
            mean[c][i]=tmp[i]/data.size();
        }
    }
    for(int c=0;c<cluster;c++)
    {
        for(int i=0;i<dim;i++)
        {
            mean[c][i]=rand()*10;
        }
    }*/
    
}
void GMM::calculate()
{
    int cc=1;
    float nowDiff=this->diff();
    log.push_back("開始計算");
    //calBeta();
    while(nowDiff>convergence)
    {
        time(&now);
        log.push_back("開始第"+to_string (cc)+"次更新於 "+(string)ctime(&now));
        printLog(log);
        this->calBeta();
        printLog(log);
        
        this->updateMean();
        printLog(log);
        this->updateVar();
        printLog(log);
        this->updatePart();
        
        
        nowDiff=this->diff();
        cout<<"差值為："<<nowDiff<<endl;
        log.push_back("差值為："+to_string(nowDiff));
        
        log.push_back("############################");
        printLog(log);
        
        cc++;
        
    }
    time(&now);
    log.push_back("計算完畢於 "+(string)ctime(&now));
}
float GMM::diff()
{
    //cout<<"計算"<<endl;
    static int counter=0;
    float d=0;
    vector<float> now;
    now.resize(0);
    if(counter>0)
    {
        //載入mean
        for(int c=0;c<cluster;c++)
        {
            for(int i=0;i<dim;i++)
            {
                now.push_back(mean[c][i]);
            }
        }
        //載入var
        for(int c=0;c<cluster;c++)
        {
            now.push_back(var[c]);
        }
        //載入part
        for(int c=0;c<cluster;c++)
        {
            now.push_back(part[c]);
        }
        prePara.resize(now.size());
        if(now.size()!=prePara.size())
            cout<<"ERROR"<<endl;
        //計算差值
        for(int i=0;i<now.size();i++)
        {
            d+=abs(now[i]-prePara[i]);
        }
        
    }
    prePara.resize(now.size());
    for(int i=0;i<now.size();i++)
    {
        prePara[i]=now[i];
    }
    
    if(counter==0)
    {
        counter++;
        return 1000000;
    }
    counter++;
    cout<<counter<<"\t";
    return d;
}
void GMM::calBeta()
{
    time(&now);
    log.push_back("開始更新Beta參數於 "+(string)ctime(&now));
    for(int i=0;i<data.size();i++)
    {
        for(int c=0;c<cluster;c++)
        {
            nowBeta[i][c]=beta(i,c);
        }
    }
    time(&now);
    log.push_back("結束更新Beta參數於 "+(string)ctime(&now));
}
void GMM::updateMean()
{
    time(&now);
    log.push_back("開始更新群心於 "+(string)ctime(&now));
    int l=this->data.size();
    //更新每一群
    for(int c=0;c<this->cluster;c++)
    {
        //初始化參數，分母
        float total=0;
        for(int i =0;i<l;i++)
        {
            total+=nowBeta[i][c];
        }
        vector<float>v;
        v.resize(this->dim);
        //分子
        for(int i =0;i<l;i++)
        {
            vector<float> x;
            x.resize(dim);
            //x
            for(int j =0;j<this->dim;j++)
            {
                x[j]=data[i][j];
            }
            multipleVector(x,nowBeta[i][c],this->dim);
            addVector(v ,x, this->dim);
        }
        divVector(v,total,this->dim);
        //更新mean
        for(int i=0;i<this->dim;i++)
        {
            mean[c][i]=v[i];
        }
    }
    time(&now);
    log.push_back("結束更新群心於 "+(string)ctime(&now));
}
void GMM::updateVar()
{
    time(&now);
    log.push_back("開始更新變異數於 "+(string)ctime(&now));
    int l=this->data.size();
    //更新每一群
    for(int c=0;c<cluster;c++)
    {
        //初始化參數，分母
        float total=0;
        for(int i =0;i<l;i++)
        {
            total+=nowBeta[i][c];
        }
        total*=this->dim;
        //分子
        float a=0;
        for(int i=0;i<l;i++)
        {
            float b=0;
            for(int j=0;j<dim;j++)
            {
                b+=pow(data[i][j]-mean[c][j],2);
            }
            b*=nowBeta[i][c];
            a+=b;
        }
        var[c]=a/total;
    }
    time(&now);
    log.push_back("結束更新變異數於 "+(string)ctime(&now));
}
void GMM::updatePart()
{
    time(&now);
    log.push_back("開始更新比重於 "+(string)ctime(&now));
    int l=this->data.size();

    for(int c=0;c<cluster;c++)
    {
        float res=0;
        for(int i=0;i<l;i++)
        {
            //cout<<beta(i,c)<<endl;
            res+=nowBeta[i][c];
        }
        //cout<<res;
        part[c]=(float)res/l;
    }
    time(&now);
    log.push_back("結束更新比重於 "+(string)ctime(&now));
}
float GMM::density(int f,int c)
{
    float res=0;
    vector<float> v;
    v.resize(dim);
    for(int i=0;i<dim;i++)
    {
        v[i]=data[f][i]-mean[c][i];
    }
    for(int i=0;i<dim;i++)
    {
        res+= pow(v[i],2);
    }
    res=res*(-1);
    res=res/(2*var[c]);
    res=exp(res);
    res*=pow(2*M_PI,(-1)*dim/2)*pow(var[c],(-1)*dim/2);
    return res;
}
float GMM::beta(int f,int c)
{
    float res=0;
    float total=0;
    for(int i=0;i<this->cluster;i++)
    {
        total+=this->density(f,i);
    }
    
    res=this->density(f,c)/total;
    return res;
    
}
void GMM::printPar()
{
    cout<<"################################################"<<endl;
    cout<<"各群的平均位置："<<endl;
    for(int c=0;c<cluster;c++)
    {
        for(int i=0;i<dim;i++)
        {
            cout<<mean[c][i]<<"\t";
        }
        cout<<endl;
    }
    
    cout<<"各群的變異數："<<endl;
    for(int c=0;c<cluster;c++)
    {
        cout<<var[c]<<"\t";
    }
    cout<<endl;
    
    cout<<"各群比例："<<endl;
    for(int c=0;c<cluster;c++)
    {
        cout<<part[c]<<"\t";
    }
    cout<<endl;
    cout<<"################################################"<<endl;
}
void GMM::Kmean(float std)
{
    /*
    //表示各群的筆數
    vector<int> num;
    num.resize(cluster);
    for(int i=0;i<cluster;i++)
    {
        num[i]=1;
    }*/
    //取得所有資料的平均
    vector<float> tmp;
    tmp.resize(dim);
    for(int i=0;i<data.size();i++)
    {
        for(int j=0;j<dim;j++)
        {
            tmp[j]+=data[i][j];
        }
    }
    for(int i=0;i<dim;i++)
    {
        tmp[i]/=data.size();
    }
    //依照平均值做常態分布，初始化各群群心
    time_t seed;
    time(&seed);
    srand(seed);
    for(int c=0;c<cluster;c++)
    {
        for(int i=0;i<dim;i++)
        {
            //mean[c][i]=tmp[i];
            //int tmp=rand()%1;
            mean[c][i]=(float)rand()/RAND_MAX;
            //gauss_distribution(tmp[i],tmp[i]*1);
        }
    }
    //計算變異數
    float tmpV=0;
    //tmpV.resize(dim);
    for(int i=0;i<data.size();i++)
    {
        tmpV+=distance(i,tmp);
    }
    tmpV/=data.size();
    for(int c=0;c<cluster;c++)
    {
        this->var[c]=tmpV/2;
    }
    
    
    //開始kmean分群
    
}
float GMM::distance(int x,vector<float> m)
{
    float res=0;
    for(int i=0;i<dim;i++)
    {
        res+=pow(data[x][i]-m[i],2);
    }
    //res=pow(res,1/2);
    return res;
}
float GMM::distance(int x,int c)
{
    float res=0;
    for(int i=0;i<dim;i++)
    {
        res+=pow(data[x][i]-mean[c][i],2);
    }
    //res=pow(res,1/2);
    return res;
}
void GMM::normalize()
{
    for(int d=0;d<dim;d++)
    {
        //找出極值
        float max=-100000;
        float min=100000;
        for(int i=0;i<data.size();i++)
        {
            if(data[i][d]>max)
                max=data[i][d];
            else if(data[i][d]<min)
                min=data[i][d];
        }
        float interval=max-min;
        
        for(int i=0;i<data.size();i++)
        {
            data[i][d]=(data[i][d]-min)/interval;
        }
    }
}
void GMM::save()
{
    string head;
    head+="設定群數為"+to_string(cluster)+"\n";
    head+="總樣本個數為"+to_string(data.size())+"\n";
    head+="樣本維度為"+to_string(dim)+"\n";
    head+="收斂條件為差值小於"+to_string(convergence)+"\n";
    head+="=============================\n";
    
    ofstream file("/record/save_"+to_string(cluster));
    file<<head;
    //輸出part
    string s="";
    for(int c=0;c<cluster;c++)
    {
        s+=to_string(part[c])+" ";
    }
    file<<s<<"\n";
    file<<"=============================\n";
    //輸出var
    s="";
    for(int c=0;c<cluster;c++)
    {
        s+=to_string(var[c])+" ";
    }
    file<<s<<"\n";
    file<<"=============================\n";
    //輸出mean
    s="";
    for(int c=0;c<cluster;c++)
    {
        for(int d=0;d<dim;d++)
        {
            s+=to_string(mean[c][d])+" ";
        }
        s+=+"\n";
    }
    file<<s<<"\n";
    file<<"=============================\n";
    
    
    file.close();
}
void GMM::load()
{}


//###########################################################
float gauss_distribution(float mean, float std)
{
    float u = rand() / (float)RAND_MAX;
    float v = rand() / (float)RAND_MAX;
    float x = sqrt(-2 * log(u)) * cos(2 * M_PI * v) * std + mean;
    
    return x;
}
void multipleVector(vector<float>& v, float m,int dim)
{
    for(int i=0;i<dim;i++)
    {
        v[i]*=m;
    }
}
void addVector(vector<float>& v1 ,vector<float>v2 ,int dim )
{
    for(int i=0;i<dim;i++)
    {
        v1[i]+=v2[i];
    }
}
void divVector(vector<float>& v1 ,float total,int dim )
{
    for(int i=0;i<dim;i++)
    {
        v1[i]=v1[i]/total;
    }
}
void printLog(vector<string> v){
    ofstream file("log");
    
    if (!file) return;
    
    for(int i=0;i<v.size();i++)
    {
        file<<v[i]<<"\n";
        //if(v[i]==1)
        //file<<"\n\n";
    }
    //file<<endl;
    file.close();
}
