#include "search.h"
#include "timer.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <cstring>
#include "timsort.hpp"

using namespace std;

Search::Search()
{
    n=0;
    sum=0;
    numDiffCoord=0;
    Tsum=0;
    Psum=0;
    IWhi=0;
    ITT=0;
    Ix=0;
    Iy=0;
    ILim=0;
    ITar=0;
}

void Search::start(const char *inputFile, const char *outputFile, bool enableTime)
{
    Timer clk;
    enableProcessTime=enableTime;

    ifstream inFile(inputFile);

    if(!inFile)
    {
        cerr << "Error while loading input file" << endl;
        return;
    }

    string line;
    int id;
    double x, y, r;

    getline(inFile, line);
    istringstream stream(line);
    stream >> n;

    initialize(n);

    for(int i=0; i<n; i++)
    {
        getline(inFile, line);
        istringstream dataStream(line);
        dataStream >> id >> x >> y >> r;
        corr[i]=id;
        Datx[i]=x;
        Daty[i]=y;
        DatR[i]=r;
    }

    inFile.close();

    clk.restart();

    GetNeighbourPair(n, corr, Datx, Daty, DatR, Pair);
    ttTime=clk.elapsed();


    //    fstream outFile;

    //    outFile.open(outputFile, fstream::out);

    //    int *cnt=new int[n]();
    //    memset(cnt, 0, sizeof(int)*n);

    //    for(int i=0; i<Psum; i++)
    //    {
    //        outFile << Pair[i][0] << "\t" << Pair[i][1] << endl;
    //        cnt[Pair[i][0]]++;
    //        cnt[Pair[i][1]]++;
    //    }

    //    outFile.close();

    cout << "********** Time Statistics **********" << endl;
    cout << endl;
    cout << "Discretizing x costs " << dxTime << " s" << endl;
    cout << "Discretizing y costs " << dyTime << " s" << endl;
    cout << "Constucting events costs " << ceTime << " s" << endl;
    cout << "Sorting events costs " << seTime << " s" << endl;
    cout << "Initializing arrays for scanning process costs " << iaTime << " s" << endl;
    cout << "Buidling tree costs " << btTime << " s" << endl;
    cout << "Scanning costs " << scTime << " s" << endl;
    cout << "The whole process (excluding dicretizing x) costs " << nxTime << " s" << endl;
    cout << endl;
    cout << "The searching process of " << n << " particles using single thread costs " << ttTime << " s" << endl;
    cout << endl;

    cout << "Find " << Psum << " pairs" << endl;


    //    fstream cntFile;

    //    cntFile.open("count.txt", fstream::out);

    //    for(int i=0; i<n; i++)
    //    {
    //        cntFile << i << "\t" << cnt[i] << endl;
    //    }

    //    cntFile.close();

    //    delete [] cnt;

}

void Search::initialize(int npnt)
{
    tree=new Node[2*(npnt+TREE_INCREMENT)-1]();
    event=new Event*[npnt*3]();
    tmpEvent=new Event[npnt*3]();
    Datx=new double[npnt]();
    Daty=new double[npnt]();
    DatR=new double[npnt]();
    diff=new double[npnt+1]();
    rectLeft=new int[npnt]();
    rectBottom=new int[npnt]();
    rectRight=new int[npnt]();
    rectTop=new int[npnt]();
    pointX=new int[npnt]();
    pointY=new int[npnt]();
    Pair=new NeighbourPair[npnt*PAIR_FACTOR]();
    rectBtmTime=new int[npnt];
    Fir=new int[npnt];
    Next=new int[npnt]();
    ITime=new int[npnt]();
    corr=new int[npnt]();
    FirSp=new int[npnt]();
    NextSp=new int[npnt*3]();
    NextS=new int[npnt*3]();

    double size=sizeof(int)*(18*npnt+1)+sizeof(double)*(4*npnt+2*PAIR_FACTOR*npnt)+sizeof(Node)*(2*(npnt+TREE_INCREMENT)-1)+sizeof(Event)*(npnt*3);
    cout << "**********Memory**********" << endl;
    cout << size/(1024*1024) << " MB" << endl;
    cout << endl;
}


void Search::GetNeighbourPair(int n, int corr[], double Datx[], double Daty[], double DatR[], int Pair[][2])
{
    Timer clk;

    discretize(Datx, DatR, rectLeft, pointX, rectRight, Horizontal);
    dxTime=clk.elapsed();

    Timer nxClk;
    clk.restart();
    discretize(Daty, DatR, rectBottom, pointY, rectTop, Vertical);
    dyTime=clk.elapsed();

    //#pragma omp parallel for
    clk.restart();
    Event *temp;
    for(int i=0; i<n*3; i++)
    {
        temp=&(tmpEvent[i]);
        temp->X= i%3==0 ? rectLeft[i/3] : (i%3==1 ? pointX[i/3] : rectRight[i/3]);
        temp->Type=i%3;
        temp->bottomPoint= (i%3==1) ? pointY[i/3] : rectBottom[i/3];
        temp->topPoint= (i%3==1) ? pointY[i/3] : rectTop[i/3];
        temp->Which=i/3;
    }
    ceTime=clk.elapsed();
    //printf("%lf\n",(double)(Clo2-Start)/CLOCKS_PER_SEC);

    //sort(c, c+3*n, Cm);
    clk.restart();
    bucketSort();
    seTime=clk.elapsed();

    //printf("%lf\n",(double)(Clo3-Clo2)/CLOCKS_PER_SEC);
    //#pragma omp parallel for
    clk.restart();
    //    for(int i=0; i<n; i++)
    //    {
    //        rectBtmTime[i]=-1;
    //        Fir[i]=-1;
    //    }
    memset(rectBtmTime, -1, sizeof(int)*n);
    memset(Fir, -1, sizeof(int)*n);
    //!Caution #pragma omp parallel for
    for (int i=0; i<n*3; i++)
    {
        if(rectBtmTime[event[i]->Which]==-1)
            rectBtmTime[event[i]->Which]=i;

        if(event[i]->Type%3==1)
            ITime[event[i]->Which]=i;
    }
    iaTime=clk.elapsed();

    clk.restart();
    BuildTree(1, numDiffCoord+TREE_INCREMENT);
    btTime=clk.elapsed();

    clk.restart();
    for(int i=0; i<n*3; i++)
    {
        if(event[i]->Type==1)
        {
            IWhi=event[i]->Which;
            ITT=i;
            Ix=event[i]->bottomPoint;
            int tmp;
            Insert(1, tmp);
        }
        else if(event[i]->Type==2)
        {
            ILim= (event[i]->Type==2) ? rectBtmTime[event[i]->Which] : i-1;
            Ix=event[i]->bottomPoint;
            Iy=event[i]->topPoint;
            IWhi=event[i]->Which;
            if(tree[1].Last>=ILim)
                GetPair(1);
        }
    }
    scTime=clk.elapsed();

    //  int nPair=0;
    //  for(int i=1; i<=Psum; i++)
    //  {
    //    if(Dist(Datx[Pair[i][0]]-Datx[Pair[i][1]], Daty[Pair[i][0]]-Daty[Pair[i][1]]) < (DatR[Pair[i][0]]+DatR[Pair[i][1]])/2+FLT_EPSILON)
    //    {
    //      Pair[++nPair][0] = Pair[i][0];
    //      Pair[nPair][1] = Pair[i][1];
    //    }
    //  }
    //  Psum=nPair;

    //cout << Psum << endl;
    nxTime=nxClk.elapsed();
    //printf("%lf\n",(double)(Finish-Clo3)/CLOCKS_PER_SEC);
}

void Search::discretize(double coord[], double r[], int left[], int point[], int right[], DiscreteMode Mod)
{
    Timer clk;
    diff[0]=-1e10;

    //#pragma omp parallel for
    for(int i=0; i<n; i++)
    {
        diff[i+1]=coord[i];
    }
    if(enableProcessTime)
    {
        cout << "init " << clk.elapsed() << endl;
        clk.restart();
    }

    sort(diff, diff+n+1);
    //    timsort(diff, diff+n+1, less<double>());
    if(enableProcessTime)
    {
        cout << "sort " << clk.elapsed() << endl;
        clk.restart();
    }

    numDiffCoord=0;

    for(int i=1; i<=n; i++)
    {
        if(diff[i]-diff[i-1]>FLT_EPSILON)
            diff[++numDiffCoord]=diff[i];
    }
    if(enableProcessTime)
    {
        cout << "diff " << clk.elapsed() << endl;
        clk.restart();
    }

    //#pragma omp parallel for
    for(int i=0; i<n; i++)
    {
        left[i]=Find(0, numDiffCoord, coord[i]-r[i]-FLT_EPSILON*2)+Mod;
        point[i]=Find(0, numDiffCoord, coord[i]);
        right[i]=Find(0, numDiffCoord, coord[i]+r[i]);
    }
    if(enableProcessTime)
    {
        cout << "dis " << clk.elapsed() << endl;
    }
}

int Search::Find(int min, int max, double target)
{
    int lowerBnd=min;
    int upperBnd=max;

    if(fabs(diff[(min+max)/2]-target)<FLT_EPSILON)
        lowerBnd=(min+max)/2;
    else
    {
        while(lowerBnd<upperBnd)
        {
            if (diff[(lowerBnd+upperBnd)/2+1]-target>FLT_EPSILON)
                upperBnd=(lowerBnd+upperBnd)/2;
            else
                lowerBnd=(lowerBnd+upperBnd)/2+1;
        }
    }
    return lowerBnd;
}

void Search::bucketSort()
{
    //    int *FirSp=new int[MAXN];
    //    int *NextSp=new int[MAXN*3];
    //    int *NextS=new int[MAXN*3];
    FirS[1]=-1;
    FirS[2]=-1;

    Timer clk;
    for(int i=0; i<3*n; i++)
    {
        int tmp= tmpEvent[i].Type==1 ? 1 : 2;
        NextS[i]=FirS[tmp];
        FirS[tmp]=i;
    }
    if(enableProcessTime)
    {
        cout << "sort a " << clk.elapsed() << endl;
        clk.restart();
    }

    memset(FirSp, -1, sizeof(int)*(n+1));
    //    for(int i=0; i<n+1; i++)
    //    {
    //        FirSp[i]=-1;
    //    }
    if(enableProcessTime)
    {
        cout << "sort init " << clk.elapsed() << endl;
        clk.restart();
    }

    for(int i=2; i>0; i--)
    {
        for(int j=FirS[i]; j!=-1; j=NextS[j])
        {
            int tmp=tmpEvent[j].X;
            NextSp[j]=FirSp[tmp];
            FirSp[tmp]=j;
        }
    }
    if(enableProcessTime)
    {
        cout << "sort b " << clk.elapsed() << endl;
        clk.restart();
    }

    int sum=-1;
    for(int i=0; i<=n; i++)
    {
        for(int j=FirSp[i]; j!=-1; j=NextSp[j])
        {
            event[++sum]=&(tmpEvent[j]);
        }
    }
    if(enableProcessTime)
    {
        cout << "sort c " << clk.elapsed() << endl;
    }

    //    delete [] FirSp;
    //    delete [] NextSp;
    //    delete [] NextS;
}

void Search::BuildTree(int x, int y)
{
    int now=++Tsum;
    tnt.x=x;
    tnt.y=y;
    tnt.sum=0;
    tnt.Last=-1;
    tnt.SecL=-1;
    tnt.toLeaf=-1;

    if(x==y)
        return;

    tnt.lson=Tsum+1;
    BuildTree(x, (x+y)/2);

    tnt.rson=Tsum+1;
    BuildTree((x+y)/2+1, y);
}

void Search::Insert(int now, int &target)
{
    tnt.sum++;
    tnt.SecL=tnt.Last;
    tnt.Last=ITT;

    if(tnt.x==tnt.y)
    {
        Next[IWhi]=Fir[tnt.x];
        Fir[tnt.x]=IWhi;
        tnt.toLeaf=now;
        target=now;
        return;
    }

    if(Ix<=tlson.y)
        Insert(tnt.lson, target);
    else
        Insert(tnt.rson, target);

    tnt.toLeaf=target;
}

void Search::GetPair(int now)
{
    if(tnt.x==tnt.y)
    {
        for(int tmp=Fir[tnt.x]; (tmp!=-1) && (ITime[tmp]>=ILim); tmp=Next[tmp])
        {
            //if(IWhi<tmp && (Dist(Datx[IWhi]-Datx[tmp], Daty[IWhi]-Daty[tmp]) < ((DatR[IWhi]+DatR[tmp])/2+FLT_EPSILON)))
            if(IWhi!=tmp)
            {
                Pair[Psum][0]=IWhi;
                Pair[Psum][1]=tmp;
                Psum++;
                assert(Psum<=n*PAIR_FACTOR);
            }
        }
        return;
    }

    if(tnt.SecL<ILim && tleaf.x>=Ix && tleaf.x<=Iy)
    {
        GetPair(tnt.toLeaf);
        return;
    }

    if(tlson.Last>=ILim && tlson.y>=Ix)
        GetPair(tnt.lson);

    if(trson.Last>=ILim && trson.x<=Iy)
        GetPair(tnt.rson);
}

double Search::Dist(double x,double y)
{
    return sqrt(x*x+y*y);
}

int Search::pairNum()
{
    return Psum;
}

Search::~Search()
{
    delete [] tree;
    delete [] event;
    delete [] tmpEvent;
    delete [] Datx;
    delete [] Daty;
    delete [] DatR;
    delete [] diff;
    delete [] rectLeft;
    delete [] rectBottom;
    delete [] rectRight;
    delete [] rectTop;
    delete [] pointX;
    delete [] pointY;
    delete [] Pair;
    delete [] rectBtmTime;
    delete [] Fir;
    delete [] Next;
    delete [] ITime;
    delete [] corr;
    delete [] NextS;
    delete [] FirSp;
    delete [] NextSp;

}
