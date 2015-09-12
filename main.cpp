#include "geninput.h"
#include "search.h"
#include "neighbour_search.h"
#include "direct.h"
#include "link_list_algorithm.h"

#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
    string mode=string(argv[1]);
    if(mode=="-gen")
    {
        if(argc == 5)
        {
            GenInput generator;
            string genTypeStr = string(argv[4]);
            generator.generate(atoi(argv[2]), atoi(argv[3]), genTypeStr);
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
    else if (mode == "-direct")
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
            Direct direct;
            direct.search(argv[2], argv[3], enableTime);
        }
        else
        {
            cout << "Error while generating file: please check the arguments!" << endl;
            return -1;
        }
    }
    else if (mode == "-linklist")
    {
        bool enableTime = false;
        if (argc >= 3)
        {
            if(argc == 4)
            {
                string enableTimeStr=string(argv[5]);
                if(enableTimeStr=="--with-process-time")
                    enableTime=true;
            }
            LinkListAlgorithm lla;
            lla.Search(argv[2]);
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
