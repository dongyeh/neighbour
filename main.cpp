#include "geninput.h"
#include "search.h"
#include "neighbour_search.h"

#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
    string mode=string(argv[1]);
    if(mode=="-gen")
    {
        if(argc==4)
        {
            GenInput generator;
            string genTypeStr = string(argv[3]);
            bool genRandom = false;
            if (genTypeStr == "random")
            {
              genRandom = true;
            }
            generator.generate(atoi(argv[2]), genRandom);
        }
        else
        {
            cout << "Error while generating file: please check the arguments!" << endl;
            return -1;
        }
    }
    else if(mode=="-run")
    {
        bool enableTime=false;
        if(argc>=5)
        {
            if(argc==6)
            {
                string enableTimeStr=string(argv[5]);
                if(enableTimeStr=="--with-process-time")
                    enableTime=true;
            }

            NeighbourSearch search(atoi(argv[4]));
            search.Search(argv[2], argv[3], enableTime);
        }
        else
        {
            cout << "Error while generating file: please check the arguments!" << endl;
            return -1;
        }
    }
    else if (mode == "-one")
    {
        bool enableTime = false;
        if (argc >= 4)
        {
            if(argc == 5)
            {
                string enableTimeStr=string(argv[5]);
                if(enableTimeStr=="--with-process-time")
                    enableTime=true;
            }
            Search search;
            search.start(argv[2], argv[3], enableTime);
        }
        else
        {
            cout << "Error while generating file: please check the arguments!" << endl;
            return -1;
        }
    }
    else
    {
        cout << "Error: illegal mode - " << mode << endl;
        return -2;
    }

    return 0;

}
