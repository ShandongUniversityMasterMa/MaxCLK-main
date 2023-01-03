#include<ctime>
#include<fstream>
#include<vector>
#include<string>
#include <algorithm>
#include<iostream>
#include<map>
using namespace std;
#define S atoi(argv[2])

#define get_genepairs(gene_pairs) { \
    while(getline(ifs1, temp)) { \
        split(temp, sto,  Line); \
        if(atoi(Line[0].c_str())>=S){ \
            SIZE=Line.size(); \
            for(i=7;i<SIZE-1;i++) { \
                for(j=i+1;j<SIZE;j++) { \
                    gene_pairs[make_pair(Line[i],Line[j])]++; \
                    if(!map1[make_pair(Line[i],Line[j])]) { \
                        all_pair.push_back(make_pair(Line[i],Line[j])); \
                        map1[make_pair(Line[i],Line[j])]=1; \
                    } \
                } \
            } \
        } \
        Line.clear(); \
    } \
    ifs1.close(); \
}

inline static void split(std::string src, std::string token, vector<std::string>& vect)
{
    int nend=0;
    int nbegin=0;
    while(nend != -1)
    {
        nend = src.find_first_of(token, nbegin);
        if(nend == -1)
            vect.push_back(src.substr(nbegin, src.length()-nbegin));
        else
            vect.push_back(src.substr(nbegin, nend-nbegin));
        nbegin = nend + 1;
    }
}

int Value(int i)
{
	int j;
	if(i>0)
	{
		j=1;
	}
	else
	{
		j=0;
	}
	return j;
}
/*
  ./get_genepairs ./CovEx/ 1 ./hi0/CSmodule_2.txt ./hi0/CSmodule_3.txt ./hi0/CSmodule_4.txt ./hi0/CSmodule_5.txt ./iref0/CSmodule_2.txt ./iref0/CSmodule_3.txt  ./iref0/CSmodule_4.txt ./iref0/CSmodule_5.txt ./mult0/CSmodule_2.txt ./mult0/CSmodule_3.txt ./mult0/CSmodule_4.txt ./mult0/CSmodule_5.txt ./hi1/CSmodule_2.txt ./hi1/CSmodule_3.txt ./hi1/CSmodule_4.txt ./hi1/CSmodule_5.txt ./iref1/CSmodule_2.txt ./iref1/CSmodule_3.txt  ./iref1/CSmodule_4.txt ./iref1/CSmodule_5.txt ./mult1/CSmodule_2.txt ./mult1/CSmodule_3.txt ./mult1/CSmodule_4.txt ./mult1/CSmodule_5.txt
*/

/*
./a.out ./Out 1 \
./pan_cancer_A/hi0_CSmodule_2.txt ./pan_cancer_A/hi0_CSmodule_3.txt ./pan_cancer_A/hi0_CSmodule_4.txt ./pan_cancer_A/hi0_CSmodule_5.txt \
./pan_cancer_A/iref0_CSmodule_2.txt ./pan_cancer_A/iref0_CSmodule_3.txt ./pan_cancer_A/iref0_CSmodule_4.txt ./pan_cancer_A/iref0_CSmodule_5.txt \
./pan_cancer_A/mult0_CSmodule_2.txt ./pan_cancer_A/mult0_CSmodule_3.txt ./pan_cancer_A/mult0_CSmodule_4.txt ./pan_cancer_A/mult0_CSmodule_5.txt \
./pan_cancer_A/hi1_CSmodule_2.txt ./pan_cancer_A/hi1_CSmodule_3.txt ./pan_cancer_A/hi1_CSmodule_4.txt ./pan_cancer_A/hi1_CSmodule_5.txt \
./pan_cancer_A/iref1_CSmodule_2.txt ./pan_cancer_A/iref1_CSmodule_3.txt ./pan_cancer_A/iref1_CSmodule_4.txt ./pan_cancer_A/iref1_CSmodule_5.txt \
./pan_cancer_A/mult1_CSmodule_2.txt ./pan_cancer_A/mult1_CSmodule_3.txt ./pan_cancer_A/mult1_CSmodule_4.txt ./pan_cancer_A/mult1_CSmodule_5.txt
*/
int main (int argc, char *argv[])
{
	int i=0,j=0,SIZE;
	int hi,iref,mult;
	vector<int> Num,int_line,int_line1;
	vector<string> Line,All_gene,T3_gene,Matrix_line,line1,line2,linkers;
	map<pair<string,string>,int> map1;
	map<pair<string, string>, int> himap_l1_u1_2, himap_l1_u1_3, himap_l1_u1_4, himap_l1_u2_2, himap_l1_u2_3, himap_l1_u2_4, himap_l1_u3_2, himap_l1_u3_3, himap_l1_u3_4, himap_l1_u4_2, himap_l1_u4_3, himap_l1_u4_4, himap_l1_u5_2, himap_l1_u5_3, himap_l1_u5_4, himap_l2_u1_2, himap_l2_u1_3, himap_l2_u1_4, himap_l2_u2_2, himap_l2_u2_3, himap_l2_u2_4, himap_l2_u3_2, himap_l2_u3_3, himap_l2_u3_4, himap_l2_u4_2, himap_l2_u4_3, himap_l2_u4_4, himap_l2_u5_2, himap_l2_u5_3, himap_l2_u5_4, himap_l3_u1_2, himap_l3_u1_3, himap_l3_u1_4, himap_l3_u2_2, himap_l3_u2_3, himap_l3_u2_4, himap_l3_u3_2, himap_l3_u3_3, himap_l3_u3_4, himap_l3_u4_2, himap_l3_u4_3, himap_l3_u4_4, himap_l3_u5_2, himap_l3_u5_3, himap_l3_u5_4, himap_l4_u1_2, himap_l4_u1_3, himap_l4_u1_4, himap_l4_u2_2, himap_l4_u2_3, himap_l4_u2_4, himap_l4_u3_2, himap_l4_u3_3, himap_l4_u3_4, himap_l4_u4_2, himap_l4_u4_3, himap_l4_u4_4, himap_l4_u5_2, himap_l4_u5_3, himap_l4_u5_4, himap_l5_u1_2, himap_l5_u1_3, himap_l5_u1_4, himap_l5_u2_2, himap_l5_u2_3, himap_l5_u2_4, himap_l5_u3_2, himap_l5_u3_3, himap_l5_u3_4, himap_l5_u4_2, himap_l5_u4_3, himap_l5_u4_4, himap_l5_u5_2, himap_l5_u5_3, himap_l5_u5_4, himap_l6_u1_2, himap_l6_u1_3, himap_l6_u1_4, himap_l6_u2_2, himap_l6_u2_3, himap_l6_u2_4, himap_l6_u3_2, himap_l6_u3_3, himap_l6_u3_4, himap_l6_u4_2, himap_l6_u4_3, himap_l6_u4_4, himap_l6_u5_2, himap_l6_u5_3, himap_l6_u5_4;
	map<pair<string, string>, int> irefmap_l1_u1_2, irefmap_l1_u1_3, irefmap_l1_u1_4, irefmap_l1_u2_2, irefmap_l1_u2_3, irefmap_l1_u2_4, irefmap_l1_u3_2, irefmap_l1_u3_3, irefmap_l1_u3_4, irefmap_l1_u4_2, irefmap_l1_u4_3, irefmap_l1_u4_4, irefmap_l1_u5_2, irefmap_l1_u5_3, irefmap_l1_u5_4, irefmap_l2_u1_2, irefmap_l2_u1_3, irefmap_l2_u1_4, irefmap_l2_u2_2, irefmap_l2_u2_3, irefmap_l2_u2_4, irefmap_l2_u3_2, irefmap_l2_u3_3, irefmap_l2_u3_4, irefmap_l2_u4_2, irefmap_l2_u4_3, irefmap_l2_u4_4, irefmap_l2_u5_2, irefmap_l2_u5_3, irefmap_l2_u5_4, irefmap_l3_u1_2, irefmap_l3_u1_3, irefmap_l3_u1_4, irefmap_l3_u2_2, irefmap_l3_u2_3, irefmap_l3_u2_4, irefmap_l3_u3_2, irefmap_l3_u3_3, irefmap_l3_u3_4, irefmap_l3_u4_2, irefmap_l3_u4_3, irefmap_l3_u4_4, irefmap_l3_u5_2, irefmap_l3_u5_3, irefmap_l3_u5_4, irefmap_l4_u1_2, irefmap_l4_u1_3, irefmap_l4_u1_4, irefmap_l4_u2_2, irefmap_l4_u2_3, irefmap_l4_u2_4, irefmap_l4_u3_2, irefmap_l4_u3_3, irefmap_l4_u3_4, irefmap_l4_u4_2, irefmap_l4_u4_3, irefmap_l4_u4_4, irefmap_l4_u5_2, irefmap_l4_u5_3, irefmap_l4_u5_4, irefmap_l5_u1_2, irefmap_l5_u1_3, irefmap_l5_u1_4, irefmap_l5_u2_2, irefmap_l5_u2_3, irefmap_l5_u2_4, irefmap_l5_u3_2, irefmap_l5_u3_3, irefmap_l5_u3_4, irefmap_l5_u4_2, irefmap_l5_u4_3, irefmap_l5_u4_4, irefmap_l5_u5_2, irefmap_l5_u5_3, irefmap_l5_u5_4, irefmap_l6_u1_2, irefmap_l6_u1_3, irefmap_l6_u1_4, irefmap_l6_u2_2, irefmap_l6_u2_3, irefmap_l6_u2_4, irefmap_l6_u3_2, irefmap_l6_u3_3, irefmap_l6_u3_4, irefmap_l6_u4_2, irefmap_l6_u4_3, irefmap_l6_u4_4, irefmap_l6_u5_2, irefmap_l6_u5_3, irefmap_l6_u5_4;
	map<pair<string, string>, int> multmap_l1_u1_2, multmap_l1_u1_3, multmap_l1_u1_4, multmap_l1_u2_2, multmap_l1_u2_3, multmap_l1_u2_4, multmap_l1_u3_2, multmap_l1_u3_3, multmap_l1_u3_4, multmap_l1_u4_2, multmap_l1_u4_3, multmap_l1_u4_4, multmap_l1_u5_2, multmap_l1_u5_3, multmap_l1_u5_4, multmap_l2_u1_2, multmap_l2_u1_3, multmap_l2_u1_4, multmap_l2_u2_2, multmap_l2_u2_3, multmap_l2_u2_4, multmap_l2_u3_2, multmap_l2_u3_3, multmap_l2_u3_4, multmap_l2_u4_2, multmap_l2_u4_3, multmap_l2_u4_4, multmap_l2_u5_2, multmap_l2_u5_3, multmap_l2_u5_4, multmap_l3_u1_2, multmap_l3_u1_3, multmap_l3_u1_4, multmap_l3_u2_2, multmap_l3_u2_3, multmap_l3_u2_4, multmap_l3_u3_2, multmap_l3_u3_3, multmap_l3_u3_4, multmap_l3_u4_2, multmap_l3_u4_3, multmap_l3_u4_4, multmap_l3_u5_2, multmap_l3_u5_3, multmap_l3_u5_4, multmap_l4_u1_2, multmap_l4_u1_3, multmap_l4_u1_4, multmap_l4_u2_2, multmap_l4_u2_3, multmap_l4_u2_4, multmap_l4_u3_2, multmap_l4_u3_3, multmap_l4_u3_4, multmap_l4_u4_2, multmap_l4_u4_3, multmap_l4_u4_4, multmap_l4_u5_2, multmap_l4_u5_3, multmap_l4_u5_4, multmap_l5_u1_2, multmap_l5_u1_3, multmap_l5_u1_4, multmap_l5_u2_2, multmap_l5_u2_3, multmap_l5_u2_4, multmap_l5_u3_2, multmap_l5_u3_3, multmap_l5_u3_4, multmap_l5_u4_2, multmap_l5_u4_3, multmap_l5_u4_4, multmap_l5_u5_2, multmap_l5_u5_3, multmap_l5_u5_4, multmap_l6_u1_2, multmap_l6_u1_3, multmap_l6_u1_4, multmap_l6_u2_2, multmap_l6_u2_3, multmap_l6_u2_4, multmap_l6_u3_2, multmap_l6_u3_3, multmap_l6_u3_4, multmap_l6_u4_2, multmap_l6_u4_3, multmap_l6_u4_4, multmap_l6_u5_2, multmap_l6_u5_3, multmap_l6_u5_4;
	map<pair<int,string>,int> map8,map9;
	vector<pair<string,string> > all_pair1,all_pair,pair1,pair2;
	vector<pair<int,string> > pair1_1,pair2_1,pair3;
	vector<vector<string> >  Matrix,CORE;

	string temp;
	string sto=" ";
	ofstream ofs,ofs1;
	char  szFileName1[100];
	char  szFileName2[100];
	sprintf(szFileName1,"%sMWgenepairs.txt",argv[1]); //multiple weights
	sprintf(szFileName2,"%sSWgenepairs.txt",argv[1]); //single weights
	ofs.open(szFileName1);
	ofs1.open(szFileName2);
	ifstream ifs1;

	ifs1.open(argv[3]); if (ifs1.is_open()) { get_genepairs(himap_l1_u1_2); ifs1.close(); }
	ifs1.open(argv[4]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u1_2); ifs1.close(); }
	ifs1.open(argv[5]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u1_2); ifs1.close(); }
	ifs1.open(argv[6]); if (ifs1.is_open()) { get_genepairs(himap_l1_u1_3); ifs1.close(); }
	ifs1.open(argv[7]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u1_3); ifs1.close(); }
	ifs1.open(argv[8]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u1_3); ifs1.close(); }
	ifs1.open(argv[9]); if (ifs1.is_open()) { get_genepairs(himap_l1_u1_4); ifs1.close(); }
	ifs1.open(argv[10]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u1_4); ifs1.close(); }
	ifs1.open(argv[11]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u1_4); ifs1.close(); }
	ifs1.open(argv[12]); if (ifs1.is_open()) { get_genepairs(himap_l1_u2_2); ifs1.close(); }
	ifs1.open(argv[13]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u2_2); ifs1.close(); }
	ifs1.open(argv[14]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u2_2); ifs1.close(); }
	ifs1.open(argv[15]); if (ifs1.is_open()) { get_genepairs(himap_l1_u2_3); ifs1.close(); }
	ifs1.open(argv[16]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u2_3); ifs1.close(); }
	ifs1.open(argv[17]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u2_3); ifs1.close(); }
	ifs1.open(argv[18]); if (ifs1.is_open()) { get_genepairs(himap_l1_u2_4); ifs1.close(); }
	ifs1.open(argv[19]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u2_4); ifs1.close(); }
	ifs1.open(argv[20]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u2_4); ifs1.close(); }
	ifs1.open(argv[21]); if (ifs1.is_open()) { get_genepairs(himap_l1_u3_2); ifs1.close(); }
	ifs1.open(argv[22]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u3_2); ifs1.close(); }
	ifs1.open(argv[23]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u3_2); ifs1.close(); }
	ifs1.open(argv[24]); if (ifs1.is_open()) { get_genepairs(himap_l1_u3_3); ifs1.close(); }
	ifs1.open(argv[25]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u3_3); ifs1.close(); }
	ifs1.open(argv[26]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u3_3); ifs1.close(); }
	ifs1.open(argv[27]); if (ifs1.is_open()) { get_genepairs(himap_l1_u3_4); ifs1.close(); }
	ifs1.open(argv[28]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u3_4); ifs1.close(); }
	ifs1.open(argv[29]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u3_4); ifs1.close(); }
	ifs1.open(argv[30]); if (ifs1.is_open()) { get_genepairs(himap_l1_u4_2); ifs1.close(); }
	ifs1.open(argv[31]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u4_2); ifs1.close(); }
	ifs1.open(argv[32]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u4_2); ifs1.close(); }
	ifs1.open(argv[33]); if (ifs1.is_open()) { get_genepairs(himap_l1_u4_3); ifs1.close(); }
	ifs1.open(argv[34]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u4_3); ifs1.close(); }
	ifs1.open(argv[35]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u4_3); ifs1.close(); }
	ifs1.open(argv[36]); if (ifs1.is_open()) { get_genepairs(himap_l1_u4_4); ifs1.close(); }
	ifs1.open(argv[37]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u4_4); ifs1.close(); }
	ifs1.open(argv[38]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u4_4); ifs1.close(); }
	ifs1.open(argv[39]); if (ifs1.is_open()) { get_genepairs(himap_l1_u5_2); ifs1.close(); }
	ifs1.open(argv[40]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u5_2); ifs1.close(); }
	ifs1.open(argv[41]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u5_2); ifs1.close(); }
	ifs1.open(argv[42]); if (ifs1.is_open()) { get_genepairs(himap_l1_u5_3); ifs1.close(); }
	ifs1.open(argv[43]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u5_3); ifs1.close(); }
	ifs1.open(argv[44]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u5_3); ifs1.close(); }
	ifs1.open(argv[45]); if (ifs1.is_open()) { get_genepairs(himap_l1_u5_4); ifs1.close(); }
	ifs1.open(argv[46]); if (ifs1.is_open()) { get_genepairs(irefmap_l1_u5_4); ifs1.close(); }
	ifs1.open(argv[47]); if (ifs1.is_open()) { get_genepairs(multmap_l1_u5_4); ifs1.close(); }
	ifs1.open(argv[48]); if (ifs1.is_open()) { get_genepairs(himap_l2_u1_2); ifs1.close(); }
	ifs1.open(argv[49]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u1_2); ifs1.close(); }
	ifs1.open(argv[50]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u1_2); ifs1.close(); }
	ifs1.open(argv[51]); if (ifs1.is_open()) { get_genepairs(himap_l2_u1_3); ifs1.close(); }
	ifs1.open(argv[52]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u1_3); ifs1.close(); }
	ifs1.open(argv[53]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u1_3); ifs1.close(); }
	ifs1.open(argv[54]); if (ifs1.is_open()) { get_genepairs(himap_l2_u1_4); ifs1.close(); }
	ifs1.open(argv[55]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u1_4); ifs1.close(); }
	ifs1.open(argv[56]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u1_4); ifs1.close(); }
	ifs1.open(argv[57]); if (ifs1.is_open()) { get_genepairs(himap_l2_u2_2); ifs1.close(); }
	ifs1.open(argv[58]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u2_2); ifs1.close(); }
	ifs1.open(argv[59]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u2_2); ifs1.close(); }
	ifs1.open(argv[60]); if (ifs1.is_open()) { get_genepairs(himap_l2_u2_3); ifs1.close(); }
	ifs1.open(argv[61]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u2_3); ifs1.close(); }
	ifs1.open(argv[62]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u2_3); ifs1.close(); }
	ifs1.open(argv[63]); if (ifs1.is_open()) { get_genepairs(himap_l2_u2_4); ifs1.close(); }
	ifs1.open(argv[64]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u2_4); ifs1.close(); }
	ifs1.open(argv[65]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u2_4); ifs1.close(); }
	ifs1.open(argv[66]); if (ifs1.is_open()) { get_genepairs(himap_l2_u3_2); ifs1.close(); }
	ifs1.open(argv[67]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u3_2); ifs1.close(); }
	ifs1.open(argv[68]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u3_2); ifs1.close(); }
	ifs1.open(argv[69]); if (ifs1.is_open()) { get_genepairs(himap_l2_u3_3); ifs1.close(); }
	ifs1.open(argv[70]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u3_3); ifs1.close(); }
	ifs1.open(argv[71]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u3_3); ifs1.close(); }
	ifs1.open(argv[72]); if (ifs1.is_open()) { get_genepairs(himap_l2_u3_4); ifs1.close(); }
	ifs1.open(argv[73]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u3_4); ifs1.close(); }
	ifs1.open(argv[74]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u3_4); ifs1.close(); }
	ifs1.open(argv[75]); if (ifs1.is_open()) { get_genepairs(himap_l2_u4_2); ifs1.close(); }
	ifs1.open(argv[76]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u4_2); ifs1.close(); }
	ifs1.open(argv[77]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u4_2); ifs1.close(); }
	ifs1.open(argv[78]); if (ifs1.is_open()) { get_genepairs(himap_l2_u4_3); ifs1.close(); }
	ifs1.open(argv[79]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u4_3); ifs1.close(); }
	ifs1.open(argv[80]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u4_3); ifs1.close(); }
	ifs1.open(argv[81]); if (ifs1.is_open()) { get_genepairs(himap_l2_u4_4); ifs1.close(); }
	ifs1.open(argv[82]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u4_4); ifs1.close(); }
	ifs1.open(argv[83]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u4_4); ifs1.close(); }
	ifs1.open(argv[84]); if (ifs1.is_open()) { get_genepairs(himap_l2_u5_2); ifs1.close(); }
	ifs1.open(argv[85]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u5_2); ifs1.close(); }
	ifs1.open(argv[86]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u5_2); ifs1.close(); }
	ifs1.open(argv[87]); if (ifs1.is_open()) { get_genepairs(himap_l2_u5_3); ifs1.close(); }
	ifs1.open(argv[88]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u5_3); ifs1.close(); }
	ifs1.open(argv[89]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u5_3); ifs1.close(); }
	ifs1.open(argv[90]); if (ifs1.is_open()) { get_genepairs(himap_l2_u5_4); ifs1.close(); }
	ifs1.open(argv[91]); if (ifs1.is_open()) { get_genepairs(irefmap_l2_u5_4); ifs1.close(); }
	ifs1.open(argv[92]); if (ifs1.is_open()) { get_genepairs(multmap_l2_u5_4); ifs1.close(); }
	ifs1.open(argv[93]); if (ifs1.is_open()) { get_genepairs(himap_l3_u1_2); ifs1.close(); }
	ifs1.open(argv[94]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u1_2); ifs1.close(); }
	ifs1.open(argv[95]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u1_2); ifs1.close(); }
	ifs1.open(argv[96]); if (ifs1.is_open()) { get_genepairs(himap_l3_u1_3); ifs1.close(); }
	ifs1.open(argv[97]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u1_3); ifs1.close(); }
	ifs1.open(argv[98]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u1_3); ifs1.close(); }
	ifs1.open(argv[99]); if (ifs1.is_open()) { get_genepairs(himap_l3_u1_4); ifs1.close(); }
	ifs1.open(argv[100]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u1_4); ifs1.close(); }
	ifs1.open(argv[101]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u1_4); ifs1.close(); }
	ifs1.open(argv[102]); if (ifs1.is_open()) { get_genepairs(himap_l3_u2_2); ifs1.close(); }
	ifs1.open(argv[103]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u2_2); ifs1.close(); }
	ifs1.open(argv[104]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u2_2); ifs1.close(); }
	ifs1.open(argv[105]); if (ifs1.is_open()) { get_genepairs(himap_l3_u2_3); ifs1.close(); }
	ifs1.open(argv[106]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u2_3); ifs1.close(); }
	ifs1.open(argv[107]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u2_3); ifs1.close(); }
	ifs1.open(argv[108]); if (ifs1.is_open()) { get_genepairs(himap_l3_u2_4); ifs1.close(); }
	ifs1.open(argv[109]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u2_4); ifs1.close(); }
	ifs1.open(argv[110]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u2_4); ifs1.close(); }
	ifs1.open(argv[111]); if (ifs1.is_open()) { get_genepairs(himap_l3_u3_2); ifs1.close(); }
	ifs1.open(argv[112]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u3_2); ifs1.close(); }
	ifs1.open(argv[113]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u3_2); ifs1.close(); }
	ifs1.open(argv[114]); if (ifs1.is_open()) { get_genepairs(himap_l3_u3_3); ifs1.close(); }
	ifs1.open(argv[115]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u3_3); ifs1.close(); }
	ifs1.open(argv[116]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u3_3); ifs1.close(); }
	ifs1.open(argv[117]); if (ifs1.is_open()) { get_genepairs(himap_l3_u3_4); ifs1.close(); }
	ifs1.open(argv[118]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u3_4); ifs1.close(); }
	ifs1.open(argv[119]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u3_4); ifs1.close(); }
	ifs1.open(argv[120]); if (ifs1.is_open()) { get_genepairs(himap_l3_u4_2); ifs1.close(); }
	ifs1.open(argv[121]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u4_2); ifs1.close(); }
	ifs1.open(argv[122]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u4_2); ifs1.close(); }
	ifs1.open(argv[123]); if (ifs1.is_open()) { get_genepairs(himap_l3_u4_3); ifs1.close(); }
	ifs1.open(argv[124]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u4_3); ifs1.close(); }
	ifs1.open(argv[125]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u4_3); ifs1.close(); }
	ifs1.open(argv[126]); if (ifs1.is_open()) { get_genepairs(himap_l3_u4_4); ifs1.close(); }
	ifs1.open(argv[127]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u4_4); ifs1.close(); }
	ifs1.open(argv[128]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u4_4); ifs1.close(); }
	ifs1.open(argv[129]); if (ifs1.is_open()) { get_genepairs(himap_l3_u5_2); ifs1.close(); }
	ifs1.open(argv[130]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u5_2); ifs1.close(); }
	ifs1.open(argv[131]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u5_2); ifs1.close(); }
	ifs1.open(argv[132]); if (ifs1.is_open()) { get_genepairs(himap_l3_u5_3); ifs1.close(); }
	ifs1.open(argv[133]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u5_3); ifs1.close(); }
	ifs1.open(argv[134]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u5_3); ifs1.close(); }
	ifs1.open(argv[135]); if (ifs1.is_open()) { get_genepairs(himap_l3_u5_4); ifs1.close(); }
	ifs1.open(argv[136]); if (ifs1.is_open()) { get_genepairs(irefmap_l3_u5_4); ifs1.close(); }
	ifs1.open(argv[137]); if (ifs1.is_open()) { get_genepairs(multmap_l3_u5_4); ifs1.close(); }
	ifs1.open(argv[138]); if (ifs1.is_open()) { get_genepairs(himap_l4_u1_2); ifs1.close(); }
	ifs1.open(argv[139]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u1_2); ifs1.close(); }
	ifs1.open(argv[140]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u1_2); ifs1.close(); }
	ifs1.open(argv[141]); if (ifs1.is_open()) { get_genepairs(himap_l4_u1_3); ifs1.close(); }
	ifs1.open(argv[142]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u1_3); ifs1.close(); }
	ifs1.open(argv[143]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u1_3); ifs1.close(); }
	ifs1.open(argv[144]); if (ifs1.is_open()) { get_genepairs(himap_l4_u1_4); ifs1.close(); }
	ifs1.open(argv[145]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u1_4); ifs1.close(); }
	ifs1.open(argv[146]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u1_4); ifs1.close(); }
	ifs1.open(argv[147]); if (ifs1.is_open()) { get_genepairs(himap_l4_u2_2); ifs1.close(); }
	ifs1.open(argv[148]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u2_2); ifs1.close(); }
	ifs1.open(argv[149]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u2_2); ifs1.close(); }
	ifs1.open(argv[150]); if (ifs1.is_open()) { get_genepairs(himap_l4_u2_3); ifs1.close(); }
	ifs1.open(argv[151]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u2_3); ifs1.close(); }
	ifs1.open(argv[152]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u2_3); ifs1.close(); }
	ifs1.open(argv[153]); if (ifs1.is_open()) { get_genepairs(himap_l4_u2_4); ifs1.close(); }
	ifs1.open(argv[154]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u2_4); ifs1.close(); }
	ifs1.open(argv[155]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u2_4); ifs1.close(); }
	ifs1.open(argv[156]); if (ifs1.is_open()) { get_genepairs(himap_l4_u3_2); ifs1.close(); }
	ifs1.open(argv[157]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u3_2); ifs1.close(); }
	ifs1.open(argv[158]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u3_2); ifs1.close(); }
	ifs1.open(argv[159]); if (ifs1.is_open()) { get_genepairs(himap_l4_u3_3); ifs1.close(); }
	ifs1.open(argv[160]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u3_3); ifs1.close(); }
	ifs1.open(argv[161]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u3_3); ifs1.close(); }
	ifs1.open(argv[162]); if (ifs1.is_open()) { get_genepairs(himap_l4_u3_4); ifs1.close(); }
	ifs1.open(argv[163]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u3_4); ifs1.close(); }
	ifs1.open(argv[164]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u3_4); ifs1.close(); }
	ifs1.open(argv[165]); if (ifs1.is_open()) { get_genepairs(himap_l4_u4_2); ifs1.close(); }
	ifs1.open(argv[166]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u4_2); ifs1.close(); }
	ifs1.open(argv[167]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u4_2); ifs1.close(); }
	ifs1.open(argv[168]); if (ifs1.is_open()) { get_genepairs(himap_l4_u4_3); ifs1.close(); }
	ifs1.open(argv[169]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u4_3); ifs1.close(); }
	ifs1.open(argv[170]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u4_3); ifs1.close(); }
	ifs1.open(argv[171]); if (ifs1.is_open()) { get_genepairs(himap_l4_u4_4); ifs1.close(); }
	ifs1.open(argv[172]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u4_4); ifs1.close(); }
	ifs1.open(argv[173]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u4_4); ifs1.close(); }
	ifs1.open(argv[174]); if (ifs1.is_open()) { get_genepairs(himap_l4_u5_2); ifs1.close(); }
	ifs1.open(argv[175]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u5_2); ifs1.close(); }
	ifs1.open(argv[176]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u5_2); ifs1.close(); }
	ifs1.open(argv[177]); if (ifs1.is_open()) { get_genepairs(himap_l4_u5_3); ifs1.close(); }
	ifs1.open(argv[178]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u5_3); ifs1.close(); }
	ifs1.open(argv[179]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u5_3); ifs1.close(); }
	ifs1.open(argv[180]); if (ifs1.is_open()) { get_genepairs(himap_l4_u5_4); ifs1.close(); }
	ifs1.open(argv[181]); if (ifs1.is_open()) { get_genepairs(irefmap_l4_u5_4); ifs1.close(); }
	ifs1.open(argv[182]); if (ifs1.is_open()) { get_genepairs(multmap_l4_u5_4); ifs1.close(); }
	ifs1.open(argv[183]); if (ifs1.is_open()) { get_genepairs(himap_l5_u1_2); ifs1.close(); }
	ifs1.open(argv[184]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u1_2); ifs1.close(); }
	ifs1.open(argv[185]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u1_2); ifs1.close(); }
	ifs1.open(argv[186]); if (ifs1.is_open()) { get_genepairs(himap_l5_u1_3); ifs1.close(); }
	ifs1.open(argv[187]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u1_3); ifs1.close(); }
	ifs1.open(argv[188]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u1_3); ifs1.close(); }
	ifs1.open(argv[189]); if (ifs1.is_open()) { get_genepairs(himap_l5_u1_4); ifs1.close(); }
	ifs1.open(argv[190]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u1_4); ifs1.close(); }
	ifs1.open(argv[191]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u1_4); ifs1.close(); }
	ifs1.open(argv[192]); if (ifs1.is_open()) { get_genepairs(himap_l5_u2_2); ifs1.close(); }
	ifs1.open(argv[193]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u2_2); ifs1.close(); }
	ifs1.open(argv[194]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u2_2); ifs1.close(); }
	ifs1.open(argv[195]); if (ifs1.is_open()) { get_genepairs(himap_l5_u2_3); ifs1.close(); }
	ifs1.open(argv[196]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u2_3); ifs1.close(); }
	ifs1.open(argv[197]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u2_3); ifs1.close(); }
	ifs1.open(argv[198]); if (ifs1.is_open()) { get_genepairs(himap_l5_u2_4); ifs1.close(); }
	ifs1.open(argv[199]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u2_4); ifs1.close(); }
	ifs1.open(argv[200]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u2_4); ifs1.close(); }
	ifs1.open(argv[201]); if (ifs1.is_open()) { get_genepairs(himap_l5_u3_2); ifs1.close(); }
	ifs1.open(argv[202]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u3_2); ifs1.close(); }
	ifs1.open(argv[203]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u3_2); ifs1.close(); }
	ifs1.open(argv[204]); if (ifs1.is_open()) { get_genepairs(himap_l5_u3_3); ifs1.close(); }
	ifs1.open(argv[205]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u3_3); ifs1.close(); }
	ifs1.open(argv[206]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u3_3); ifs1.close(); }
	ifs1.open(argv[207]); if (ifs1.is_open()) { get_genepairs(himap_l5_u3_4); ifs1.close(); }
	ifs1.open(argv[208]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u3_4); ifs1.close(); }
	ifs1.open(argv[209]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u3_4); ifs1.close(); }
	ifs1.open(argv[210]); if (ifs1.is_open()) { get_genepairs(himap_l5_u4_2); ifs1.close(); }
	ifs1.open(argv[211]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u4_2); ifs1.close(); }
	ifs1.open(argv[212]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u4_2); ifs1.close(); }
	ifs1.open(argv[213]); if (ifs1.is_open()) { get_genepairs(himap_l5_u4_3); ifs1.close(); }
	ifs1.open(argv[214]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u4_3); ifs1.close(); }
	ifs1.open(argv[215]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u4_3); ifs1.close(); }
	ifs1.open(argv[216]); if (ifs1.is_open()) { get_genepairs(himap_l5_u4_4); ifs1.close(); }
	ifs1.open(argv[217]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u4_4); ifs1.close(); }
	ifs1.open(argv[218]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u4_4); ifs1.close(); }
	ifs1.open(argv[219]); if (ifs1.is_open()) { get_genepairs(himap_l5_u5_2); ifs1.close(); }
	ifs1.open(argv[220]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u5_2); ifs1.close(); }
	ifs1.open(argv[221]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u5_2); ifs1.close(); }
	ifs1.open(argv[222]); if (ifs1.is_open()) { get_genepairs(himap_l5_u5_3); ifs1.close(); }
	ifs1.open(argv[223]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u5_3); ifs1.close(); }
	ifs1.open(argv[224]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u5_3); ifs1.close(); }
	ifs1.open(argv[225]); if (ifs1.is_open()) { get_genepairs(himap_l5_u5_4); ifs1.close(); }
	ifs1.open(argv[226]); if (ifs1.is_open()) { get_genepairs(irefmap_l5_u5_4); ifs1.close(); }
	ifs1.open(argv[227]); if (ifs1.is_open()) { get_genepairs(multmap_l5_u5_4); ifs1.close(); }
	ifs1.open(argv[228]); if (ifs1.is_open()) { get_genepairs(himap_l6_u1_2); ifs1.close(); }
	ifs1.open(argv[229]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u1_2); ifs1.close(); }
	ifs1.open(argv[230]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u1_2); ifs1.close(); }
	ifs1.open(argv[231]); if (ifs1.is_open()) { get_genepairs(himap_l6_u1_3); ifs1.close(); }
	ifs1.open(argv[232]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u1_3); ifs1.close(); }
	ifs1.open(argv[233]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u1_3); ifs1.close(); }
	ifs1.open(argv[234]); if (ifs1.is_open()) { get_genepairs(himap_l6_u1_4); ifs1.close(); }
	ifs1.open(argv[235]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u1_4); ifs1.close(); }
	ifs1.open(argv[236]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u1_4); ifs1.close(); }
	ifs1.open(argv[237]); if (ifs1.is_open()) { get_genepairs(himap_l6_u2_2); ifs1.close(); }
	ifs1.open(argv[238]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u2_2); ifs1.close(); }
	ifs1.open(argv[239]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u2_2); ifs1.close(); }
	ifs1.open(argv[240]); if (ifs1.is_open()) { get_genepairs(himap_l6_u2_3); ifs1.close(); }
	ifs1.open(argv[241]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u2_3); ifs1.close(); }
	ifs1.open(argv[242]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u2_3); ifs1.close(); }
	ifs1.open(argv[243]); if (ifs1.is_open()) { get_genepairs(himap_l6_u2_4); ifs1.close(); }
	ifs1.open(argv[244]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u2_4); ifs1.close(); }
	ifs1.open(argv[245]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u2_4); ifs1.close(); }
	ifs1.open(argv[246]); if (ifs1.is_open()) { get_genepairs(himap_l6_u3_2); ifs1.close(); }
	ifs1.open(argv[247]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u3_2); ifs1.close(); }
	ifs1.open(argv[248]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u3_2); ifs1.close(); }
	ifs1.open(argv[249]); if (ifs1.is_open()) { get_genepairs(himap_l6_u3_3); ifs1.close(); }
	ifs1.open(argv[250]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u3_3); ifs1.close(); }
	ifs1.open(argv[251]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u3_3); ifs1.close(); }
	ifs1.open(argv[252]); if (ifs1.is_open()) { get_genepairs(himap_l6_u3_4); ifs1.close(); }
	ifs1.open(argv[253]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u3_4); ifs1.close(); }
	ifs1.open(argv[254]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u3_4); ifs1.close(); }
	ifs1.open(argv[255]); if (ifs1.is_open()) { get_genepairs(himap_l6_u4_2); ifs1.close(); }
	ifs1.open(argv[256]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u4_2); ifs1.close(); }
	ifs1.open(argv[257]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u4_2); ifs1.close(); }
	ifs1.open(argv[258]); if (ifs1.is_open()) { get_genepairs(himap_l6_u4_3); ifs1.close(); }
	ifs1.open(argv[259]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u4_3); ifs1.close(); }
	ifs1.open(argv[260]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u4_3); ifs1.close(); }
	ifs1.open(argv[261]); if (ifs1.is_open()) { get_genepairs(himap_l6_u4_4); ifs1.close(); }
	ifs1.open(argv[262]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u4_4); ifs1.close(); }
	ifs1.open(argv[263]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u4_4); ifs1.close(); }
	ifs1.open(argv[264]); if (ifs1.is_open()) { get_genepairs(himap_l6_u5_2); ifs1.close(); }
	ifs1.open(argv[265]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u5_2); ifs1.close(); }
	ifs1.open(argv[266]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u5_2); ifs1.close(); }
	ifs1.open(argv[267]); if (ifs1.is_open()) { get_genepairs(himap_l6_u5_3); ifs1.close(); }
	ifs1.open(argv[268]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u5_3); ifs1.close(); }
	ifs1.open(argv[269]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u5_3); ifs1.close(); }
	ifs1.open(argv[270]); if (ifs1.is_open()) { get_genepairs(himap_l6_u5_4); ifs1.close(); }
	ifs1.open(argv[271]); if (ifs1.is_open()) { get_genepairs(irefmap_l6_u5_4); ifs1.close(); }
	ifs1.open(argv[272]); if (ifs1.is_open()) { get_genepairs(multmap_l6_u5_4); ifs1.close(); }

	SIZE=all_pair.size();
	for(i=0;i<SIZE;i++)
	{
		hi = himap_l1_u1_2[all_pair[i]] + himap_l1_u1_3[all_pair[i]] + himap_l1_u1_4[all_pair[i]] + himap_l1_u2_2[all_pair[i]] + himap_l1_u2_3[all_pair[i]] + himap_l1_u2_4[all_pair[i]] + himap_l1_u3_2[all_pair[i]] + himap_l1_u3_3[all_pair[i]] + himap_l1_u3_4[all_pair[i]] + himap_l1_u4_2[all_pair[i]] + himap_l1_u4_3[all_pair[i]] + himap_l1_u4_4[all_pair[i]] + himap_l1_u5_2[all_pair[i]] + himap_l1_u5_3[all_pair[i]] + himap_l1_u5_4[all_pair[i]] + himap_l2_u1_2[all_pair[i]] + himap_l2_u1_3[all_pair[i]] + himap_l2_u1_4[all_pair[i]] + himap_l2_u2_2[all_pair[i]] + himap_l2_u2_3[all_pair[i]] + himap_l2_u2_4[all_pair[i]] + himap_l2_u3_2[all_pair[i]] + himap_l2_u3_3[all_pair[i]] + himap_l2_u3_4[all_pair[i]] + himap_l2_u4_2[all_pair[i]] + himap_l2_u4_3[all_pair[i]] + himap_l2_u4_4[all_pair[i]] + himap_l2_u5_2[all_pair[i]] + himap_l2_u5_3[all_pair[i]] + himap_l2_u5_4[all_pair[i]] + himap_l3_u1_2[all_pair[i]] + himap_l3_u1_3[all_pair[i]] + himap_l3_u1_4[all_pair[i]] + himap_l3_u2_2[all_pair[i]] + himap_l3_u2_3[all_pair[i]] + himap_l3_u2_4[all_pair[i]] + himap_l3_u3_2[all_pair[i]] + himap_l3_u3_3[all_pair[i]] + himap_l3_u3_4[all_pair[i]] + himap_l3_u4_2[all_pair[i]] + himap_l3_u4_3[all_pair[i]] + himap_l3_u4_4[all_pair[i]] + himap_l3_u5_2[all_pair[i]] + himap_l3_u5_3[all_pair[i]] + himap_l3_u5_4[all_pair[i]] + himap_l4_u1_2[all_pair[i]] + himap_l4_u1_3[all_pair[i]] + himap_l4_u1_4[all_pair[i]] + himap_l4_u2_2[all_pair[i]] + himap_l4_u2_3[all_pair[i]] + himap_l4_u2_4[all_pair[i]] + himap_l4_u3_2[all_pair[i]] + himap_l4_u3_3[all_pair[i]] + himap_l4_u3_4[all_pair[i]] + himap_l4_u4_2[all_pair[i]] + himap_l4_u4_3[all_pair[i]] + himap_l4_u4_4[all_pair[i]] + himap_l4_u5_2[all_pair[i]] + himap_l4_u5_3[all_pair[i]] + himap_l4_u5_4[all_pair[i]] + himap_l5_u1_2[all_pair[i]] + himap_l5_u1_3[all_pair[i]] + himap_l5_u1_4[all_pair[i]] + himap_l5_u2_2[all_pair[i]] + himap_l5_u2_3[all_pair[i]] + himap_l5_u2_4[all_pair[i]] + himap_l5_u3_2[all_pair[i]] + himap_l5_u3_3[all_pair[i]] + himap_l5_u3_4[all_pair[i]] + himap_l5_u4_2[all_pair[i]] + himap_l5_u4_3[all_pair[i]] + himap_l5_u4_4[all_pair[i]] + himap_l5_u5_2[all_pair[i]] + himap_l5_u5_3[all_pair[i]] + himap_l5_u5_4[all_pair[i]] + himap_l6_u1_2[all_pair[i]] + himap_l6_u1_3[all_pair[i]] + himap_l6_u1_4[all_pair[i]] + himap_l6_u2_2[all_pair[i]] + himap_l6_u2_3[all_pair[i]] + himap_l6_u2_4[all_pair[i]] + himap_l6_u3_2[all_pair[i]] + himap_l6_u3_3[all_pair[i]] + himap_l6_u3_4[all_pair[i]] + himap_l6_u4_2[all_pair[i]] + himap_l6_u4_3[all_pair[i]] + himap_l6_u4_4[all_pair[i]] + himap_l6_u5_2[all_pair[i]] + himap_l6_u5_3[all_pair[i]] + himap_l6_u5_4[all_pair[i]];
		iref = irefmap_l1_u1_2[all_pair[i]] + irefmap_l1_u1_3[all_pair[i]] + irefmap_l1_u1_4[all_pair[i]] + irefmap_l1_u2_2[all_pair[i]] + irefmap_l1_u2_3[all_pair[i]] + irefmap_l1_u2_4[all_pair[i]] + irefmap_l1_u3_2[all_pair[i]] + irefmap_l1_u3_3[all_pair[i]] + irefmap_l1_u3_4[all_pair[i]] + irefmap_l1_u4_2[all_pair[i]] + irefmap_l1_u4_3[all_pair[i]] + irefmap_l1_u4_4[all_pair[i]] + irefmap_l1_u5_2[all_pair[i]] + irefmap_l1_u5_3[all_pair[i]] + irefmap_l1_u5_4[all_pair[i]] + irefmap_l2_u1_2[all_pair[i]] + irefmap_l2_u1_3[all_pair[i]] + irefmap_l2_u1_4[all_pair[i]] + irefmap_l2_u2_2[all_pair[i]] + irefmap_l2_u2_3[all_pair[i]] + irefmap_l2_u2_4[all_pair[i]] + irefmap_l2_u3_2[all_pair[i]] + irefmap_l2_u3_3[all_pair[i]] + irefmap_l2_u3_4[all_pair[i]] + irefmap_l2_u4_2[all_pair[i]] + irefmap_l2_u4_3[all_pair[i]] + irefmap_l2_u4_4[all_pair[i]] + irefmap_l2_u5_2[all_pair[i]] + irefmap_l2_u5_3[all_pair[i]] + irefmap_l2_u5_4[all_pair[i]] + irefmap_l3_u1_2[all_pair[i]] + irefmap_l3_u1_3[all_pair[i]] + irefmap_l3_u1_4[all_pair[i]] + irefmap_l3_u2_2[all_pair[i]] + irefmap_l3_u2_3[all_pair[i]] + irefmap_l3_u2_4[all_pair[i]] + irefmap_l3_u3_2[all_pair[i]] + irefmap_l3_u3_3[all_pair[i]] + irefmap_l3_u3_4[all_pair[i]] + irefmap_l3_u4_2[all_pair[i]] + irefmap_l3_u4_3[all_pair[i]] + irefmap_l3_u4_4[all_pair[i]] + irefmap_l3_u5_2[all_pair[i]] + irefmap_l3_u5_3[all_pair[i]] + irefmap_l3_u5_4[all_pair[i]] + irefmap_l4_u1_2[all_pair[i]] + irefmap_l4_u1_3[all_pair[i]] + irefmap_l4_u1_4[all_pair[i]] + irefmap_l4_u2_2[all_pair[i]] + irefmap_l4_u2_3[all_pair[i]] + irefmap_l4_u2_4[all_pair[i]] + irefmap_l4_u3_2[all_pair[i]] + irefmap_l4_u3_3[all_pair[i]] + irefmap_l4_u3_4[all_pair[i]] + irefmap_l4_u4_2[all_pair[i]] + irefmap_l4_u4_3[all_pair[i]] + irefmap_l4_u4_4[all_pair[i]] + irefmap_l4_u5_2[all_pair[i]] + irefmap_l4_u5_3[all_pair[i]] + irefmap_l4_u5_4[all_pair[i]] + irefmap_l5_u1_2[all_pair[i]] + irefmap_l5_u1_3[all_pair[i]] + irefmap_l5_u1_4[all_pair[i]] + irefmap_l5_u2_2[all_pair[i]] + irefmap_l5_u2_3[all_pair[i]] + irefmap_l5_u2_4[all_pair[i]] + irefmap_l5_u3_2[all_pair[i]] + irefmap_l5_u3_3[all_pair[i]] + irefmap_l5_u3_4[all_pair[i]] + irefmap_l5_u4_2[all_pair[i]] + irefmap_l5_u4_3[all_pair[i]] + irefmap_l5_u4_4[all_pair[i]] + irefmap_l5_u5_2[all_pair[i]] + irefmap_l5_u5_3[all_pair[i]] + irefmap_l5_u5_4[all_pair[i]] + irefmap_l6_u1_2[all_pair[i]] + irefmap_l6_u1_3[all_pair[i]] + irefmap_l6_u1_4[all_pair[i]] + irefmap_l6_u2_2[all_pair[i]] + irefmap_l6_u2_3[all_pair[i]] + irefmap_l6_u2_4[all_pair[i]] + irefmap_l6_u3_2[all_pair[i]] + irefmap_l6_u3_3[all_pair[i]] + irefmap_l6_u3_4[all_pair[i]] + irefmap_l6_u4_2[all_pair[i]] + irefmap_l6_u4_3[all_pair[i]] + irefmap_l6_u4_4[all_pair[i]] + irefmap_l6_u5_2[all_pair[i]] + irefmap_l6_u5_3[all_pair[i]] + irefmap_l6_u5_4[all_pair[i]];
		mult = multmap_l1_u1_2[all_pair[i]] + multmap_l1_u1_3[all_pair[i]] + multmap_l1_u1_4[all_pair[i]] + multmap_l1_u2_2[all_pair[i]] + multmap_l1_u2_3[all_pair[i]] + multmap_l1_u2_4[all_pair[i]] + multmap_l1_u3_2[all_pair[i]] + multmap_l1_u3_3[all_pair[i]] + multmap_l1_u3_4[all_pair[i]] + multmap_l1_u4_2[all_pair[i]] + multmap_l1_u4_3[all_pair[i]] + multmap_l1_u4_4[all_pair[i]] + multmap_l1_u5_2[all_pair[i]] + multmap_l1_u5_3[all_pair[i]] + multmap_l1_u5_4[all_pair[i]] + multmap_l2_u1_2[all_pair[i]] + multmap_l2_u1_3[all_pair[i]] + multmap_l2_u1_4[all_pair[i]] + multmap_l2_u2_2[all_pair[i]] + multmap_l2_u2_3[all_pair[i]] + multmap_l2_u2_4[all_pair[i]] + multmap_l2_u3_2[all_pair[i]] + multmap_l2_u3_3[all_pair[i]] + multmap_l2_u3_4[all_pair[i]] + multmap_l2_u4_2[all_pair[i]] + multmap_l2_u4_3[all_pair[i]] + multmap_l2_u4_4[all_pair[i]] + multmap_l2_u5_2[all_pair[i]] + multmap_l2_u5_3[all_pair[i]] + multmap_l2_u5_4[all_pair[i]] + multmap_l3_u1_2[all_pair[i]] + multmap_l3_u1_3[all_pair[i]] + multmap_l3_u1_4[all_pair[i]] + multmap_l3_u2_2[all_pair[i]] + multmap_l3_u2_3[all_pair[i]] + multmap_l3_u2_4[all_pair[i]] + multmap_l3_u3_2[all_pair[i]] + multmap_l3_u3_3[all_pair[i]] + multmap_l3_u3_4[all_pair[i]] + multmap_l3_u4_2[all_pair[i]] + multmap_l3_u4_3[all_pair[i]] + multmap_l3_u4_4[all_pair[i]] + multmap_l3_u5_2[all_pair[i]] + multmap_l3_u5_3[all_pair[i]] + multmap_l3_u5_4[all_pair[i]] + multmap_l4_u1_2[all_pair[i]] + multmap_l4_u1_3[all_pair[i]] + multmap_l4_u1_4[all_pair[i]] + multmap_l4_u2_2[all_pair[i]] + multmap_l4_u2_3[all_pair[i]] + multmap_l4_u2_4[all_pair[i]] + multmap_l4_u3_2[all_pair[i]] + multmap_l4_u3_3[all_pair[i]] + multmap_l4_u3_4[all_pair[i]] + multmap_l4_u4_2[all_pair[i]] + multmap_l4_u4_3[all_pair[i]] + multmap_l4_u4_4[all_pair[i]] + multmap_l4_u5_2[all_pair[i]] + multmap_l4_u5_3[all_pair[i]] + multmap_l4_u5_4[all_pair[i]] + multmap_l5_u1_2[all_pair[i]] + multmap_l5_u1_3[all_pair[i]] + multmap_l5_u1_4[all_pair[i]] + multmap_l5_u2_2[all_pair[i]] + multmap_l5_u2_3[all_pair[i]] + multmap_l5_u2_4[all_pair[i]] + multmap_l5_u3_2[all_pair[i]] + multmap_l5_u3_3[all_pair[i]] + multmap_l5_u3_4[all_pair[i]] + multmap_l5_u4_2[all_pair[i]] + multmap_l5_u4_3[all_pair[i]] + multmap_l5_u4_4[all_pair[i]] + multmap_l5_u5_2[all_pair[i]] + multmap_l5_u5_3[all_pair[i]] + multmap_l5_u5_4[all_pair[i]] + multmap_l6_u1_2[all_pair[i]] + multmap_l6_u1_3[all_pair[i]] + multmap_l6_u1_4[all_pair[i]] + multmap_l6_u2_2[all_pair[i]] + multmap_l6_u2_3[all_pair[i]] + multmap_l6_u2_4[all_pair[i]] + multmap_l6_u3_2[all_pair[i]] + multmap_l6_u3_3[all_pair[i]] + multmap_l6_u3_4[all_pair[i]] + multmap_l6_u4_2[all_pair[i]] + multmap_l6_u4_3[all_pair[i]] + multmap_l6_u4_4[all_pair[i]] + multmap_l6_u5_2[all_pair[i]] + multmap_l6_u5_3[all_pair[i]] + multmap_l6_u5_4[all_pair[i]];
		
		ofs << hi + iref + mult << " " << hi << " " << iref << " " << mult << " " ;
		ofs<<all_pair[i].first<<" "<<all_pair[i].second<<endl;
	}
	for(i=0;i<SIZE;i++)
	{
		hi = Value(himap_l1_u1_2[all_pair[i]]) + Value(himap_l1_u1_3[all_pair[i]]) + Value(himap_l1_u1_4[all_pair[i]]) + Value(himap_l1_u2_2[all_pair[i]]) + Value(himap_l1_u2_3[all_pair[i]]) + Value(himap_l1_u2_4[all_pair[i]]) + Value(himap_l1_u3_2[all_pair[i]]) + Value(himap_l1_u3_3[all_pair[i]]) + Value(himap_l1_u3_4[all_pair[i]]) + Value(himap_l1_u4_2[all_pair[i]]) + Value(himap_l1_u4_3[all_pair[i]]) + Value(himap_l1_u4_4[all_pair[i]]) + Value(himap_l1_u5_2[all_pair[i]]) + Value(himap_l1_u5_3[all_pair[i]]) + Value(himap_l1_u5_4[all_pair[i]]) + Value(himap_l2_u1_2[all_pair[i]]) + Value(himap_l2_u1_3[all_pair[i]]) + Value(himap_l2_u1_4[all_pair[i]]) + Value(himap_l2_u2_2[all_pair[i]]) + Value(himap_l2_u2_3[all_pair[i]]) + Value(himap_l2_u2_4[all_pair[i]]) + Value(himap_l2_u3_2[all_pair[i]]) + Value(himap_l2_u3_3[all_pair[i]]) + Value(himap_l2_u3_4[all_pair[i]]) + Value(himap_l2_u4_2[all_pair[i]]) + Value(himap_l2_u4_3[all_pair[i]]) + Value(himap_l2_u4_4[all_pair[i]]) + Value(himap_l2_u5_2[all_pair[i]]) + Value(himap_l2_u5_3[all_pair[i]]) + Value(himap_l2_u5_4[all_pair[i]]) + Value(himap_l3_u1_2[all_pair[i]]) + Value(himap_l3_u1_3[all_pair[i]]) + Value(himap_l3_u1_4[all_pair[i]]) + Value(himap_l3_u2_2[all_pair[i]]) + Value(himap_l3_u2_3[all_pair[i]]) + Value(himap_l3_u2_4[all_pair[i]]) + Value(himap_l3_u3_2[all_pair[i]]) + Value(himap_l3_u3_3[all_pair[i]]) + Value(himap_l3_u3_4[all_pair[i]]) + Value(himap_l3_u4_2[all_pair[i]]) + Value(himap_l3_u4_3[all_pair[i]]) + Value(himap_l3_u4_4[all_pair[i]]) + Value(himap_l3_u5_2[all_pair[i]]) + Value(himap_l3_u5_3[all_pair[i]]) + Value(himap_l3_u5_4[all_pair[i]]) + Value(himap_l4_u1_2[all_pair[i]]) + Value(himap_l4_u1_3[all_pair[i]]) + Value(himap_l4_u1_4[all_pair[i]]) + Value(himap_l4_u2_2[all_pair[i]]) + Value(himap_l4_u2_3[all_pair[i]]) + Value(himap_l4_u2_4[all_pair[i]]) + Value(himap_l4_u3_2[all_pair[i]]) + Value(himap_l4_u3_3[all_pair[i]]) + Value(himap_l4_u3_4[all_pair[i]]) + Value(himap_l4_u4_2[all_pair[i]]) + Value(himap_l4_u4_3[all_pair[i]]) + Value(himap_l4_u4_4[all_pair[i]]) + Value(himap_l4_u5_2[all_pair[i]]) + Value(himap_l4_u5_3[all_pair[i]]) + Value(himap_l4_u5_4[all_pair[i]]) + Value(himap_l5_u1_2[all_pair[i]]) + Value(himap_l5_u1_3[all_pair[i]]) + Value(himap_l5_u1_4[all_pair[i]]) + Value(himap_l5_u2_2[all_pair[i]]) + Value(himap_l5_u2_3[all_pair[i]]) + Value(himap_l5_u2_4[all_pair[i]]) + Value(himap_l5_u3_2[all_pair[i]]) + Value(himap_l5_u3_3[all_pair[i]]) + Value(himap_l5_u3_4[all_pair[i]]) + Value(himap_l5_u4_2[all_pair[i]]) + Value(himap_l5_u4_3[all_pair[i]]) + Value(himap_l5_u4_4[all_pair[i]]) + Value(himap_l5_u5_2[all_pair[i]]) + Value(himap_l5_u5_3[all_pair[i]]) + Value(himap_l5_u5_4[all_pair[i]]) + Value(himap_l6_u1_2[all_pair[i]]) + Value(himap_l6_u1_3[all_pair[i]]) + Value(himap_l6_u1_4[all_pair[i]]) + Value(himap_l6_u2_2[all_pair[i]]) + Value(himap_l6_u2_3[all_pair[i]]) + Value(himap_l6_u2_4[all_pair[i]]) + Value(himap_l6_u3_2[all_pair[i]]) + Value(himap_l6_u3_3[all_pair[i]]) + Value(himap_l6_u3_4[all_pair[i]]) + Value(himap_l6_u4_2[all_pair[i]]) + Value(himap_l6_u4_3[all_pair[i]]) + Value(himap_l6_u4_4[all_pair[i]]) + Value(himap_l6_u5_2[all_pair[i]]) + Value(himap_l6_u5_3[all_pair[i]]) + Value(himap_l6_u5_4[all_pair[i]]);
		iref = Value(irefmap_l1_u1_2[all_pair[i]]) + Value(irefmap_l1_u1_3[all_pair[i]]) + Value(irefmap_l1_u1_4[all_pair[i]]) + Value(irefmap_l1_u2_2[all_pair[i]]) + Value(irefmap_l1_u2_3[all_pair[i]]) + Value(irefmap_l1_u2_4[all_pair[i]]) + Value(irefmap_l1_u3_2[all_pair[i]]) + Value(irefmap_l1_u3_3[all_pair[i]]) + Value(irefmap_l1_u3_4[all_pair[i]]) + Value(irefmap_l1_u4_2[all_pair[i]]) + Value(irefmap_l1_u4_3[all_pair[i]]) + Value(irefmap_l1_u4_4[all_pair[i]]) + Value(irefmap_l1_u5_2[all_pair[i]]) + Value(irefmap_l1_u5_3[all_pair[i]]) + Value(irefmap_l1_u5_4[all_pair[i]]) + Value(irefmap_l2_u1_2[all_pair[i]]) + Value(irefmap_l2_u1_3[all_pair[i]]) + Value(irefmap_l2_u1_4[all_pair[i]]) + Value(irefmap_l2_u2_2[all_pair[i]]) + Value(irefmap_l2_u2_3[all_pair[i]]) + Value(irefmap_l2_u2_4[all_pair[i]]) + Value(irefmap_l2_u3_2[all_pair[i]]) + Value(irefmap_l2_u3_3[all_pair[i]]) + Value(irefmap_l2_u3_4[all_pair[i]]) + Value(irefmap_l2_u4_2[all_pair[i]]) + Value(irefmap_l2_u4_3[all_pair[i]]) + Value(irefmap_l2_u4_4[all_pair[i]]) + Value(irefmap_l2_u5_2[all_pair[i]]) + Value(irefmap_l2_u5_3[all_pair[i]]) + Value(irefmap_l2_u5_4[all_pair[i]]) + Value(irefmap_l3_u1_2[all_pair[i]]) + Value(irefmap_l3_u1_3[all_pair[i]]) + Value(irefmap_l3_u1_4[all_pair[i]]) + Value(irefmap_l3_u2_2[all_pair[i]]) + Value(irefmap_l3_u2_3[all_pair[i]]) + Value(irefmap_l3_u2_4[all_pair[i]]) + Value(irefmap_l3_u3_2[all_pair[i]]) + Value(irefmap_l3_u3_3[all_pair[i]]) + Value(irefmap_l3_u3_4[all_pair[i]]) + Value(irefmap_l3_u4_2[all_pair[i]]) + Value(irefmap_l3_u4_3[all_pair[i]]) + Value(irefmap_l3_u4_4[all_pair[i]]) + Value(irefmap_l3_u5_2[all_pair[i]]) + Value(irefmap_l3_u5_3[all_pair[i]]) + Value(irefmap_l3_u5_4[all_pair[i]]) + Value(irefmap_l4_u1_2[all_pair[i]]) + Value(irefmap_l4_u1_3[all_pair[i]]) + Value(irefmap_l4_u1_4[all_pair[i]]) + Value(irefmap_l4_u2_2[all_pair[i]]) + Value(irefmap_l4_u2_3[all_pair[i]]) + Value(irefmap_l4_u2_4[all_pair[i]]) + Value(irefmap_l4_u3_2[all_pair[i]]) + Value(irefmap_l4_u3_3[all_pair[i]]) + Value(irefmap_l4_u3_4[all_pair[i]]) + Value(irefmap_l4_u4_2[all_pair[i]]) + Value(irefmap_l4_u4_3[all_pair[i]]) + Value(irefmap_l4_u4_4[all_pair[i]]) + Value(irefmap_l4_u5_2[all_pair[i]]) + Value(irefmap_l4_u5_3[all_pair[i]]) + Value(irefmap_l4_u5_4[all_pair[i]]) + Value(irefmap_l5_u1_2[all_pair[i]]) + Value(irefmap_l5_u1_3[all_pair[i]]) + Value(irefmap_l5_u1_4[all_pair[i]]) + Value(irefmap_l5_u2_2[all_pair[i]]) + Value(irefmap_l5_u2_3[all_pair[i]]) + Value(irefmap_l5_u2_4[all_pair[i]]) + Value(irefmap_l5_u3_2[all_pair[i]]) + Value(irefmap_l5_u3_3[all_pair[i]]) + Value(irefmap_l5_u3_4[all_pair[i]]) + Value(irefmap_l5_u4_2[all_pair[i]]) + Value(irefmap_l5_u4_3[all_pair[i]]) + Value(irefmap_l5_u4_4[all_pair[i]]) + Value(irefmap_l5_u5_2[all_pair[i]]) + Value(irefmap_l5_u5_3[all_pair[i]]) + Value(irefmap_l5_u5_4[all_pair[i]]) + Value(irefmap_l6_u1_2[all_pair[i]]) + Value(irefmap_l6_u1_3[all_pair[i]]) + Value(irefmap_l6_u1_4[all_pair[i]]) + Value(irefmap_l6_u2_2[all_pair[i]]) + Value(irefmap_l6_u2_3[all_pair[i]]) + Value(irefmap_l6_u2_4[all_pair[i]]) + Value(irefmap_l6_u3_2[all_pair[i]]) + Value(irefmap_l6_u3_3[all_pair[i]]) + Value(irefmap_l6_u3_4[all_pair[i]]) + Value(irefmap_l6_u4_2[all_pair[i]]) + Value(irefmap_l6_u4_3[all_pair[i]]) + Value(irefmap_l6_u4_4[all_pair[i]]) + Value(irefmap_l6_u5_2[all_pair[i]]) + Value(irefmap_l6_u5_3[all_pair[i]]) + Value(irefmap_l6_u5_4[all_pair[i]]);
		mult = Value(multmap_l1_u1_2[all_pair[i]]) + Value(multmap_l1_u1_3[all_pair[i]]) + Value(multmap_l1_u1_4[all_pair[i]]) + Value(multmap_l1_u2_2[all_pair[i]]) + Value(multmap_l1_u2_3[all_pair[i]]) + Value(multmap_l1_u2_4[all_pair[i]]) + Value(multmap_l1_u3_2[all_pair[i]]) + Value(multmap_l1_u3_3[all_pair[i]]) + Value(multmap_l1_u3_4[all_pair[i]]) + Value(multmap_l1_u4_2[all_pair[i]]) + Value(multmap_l1_u4_3[all_pair[i]]) + Value(multmap_l1_u4_4[all_pair[i]]) + Value(multmap_l1_u5_2[all_pair[i]]) + Value(multmap_l1_u5_3[all_pair[i]]) + Value(multmap_l1_u5_4[all_pair[i]]) + Value(multmap_l2_u1_2[all_pair[i]]) + Value(multmap_l2_u1_3[all_pair[i]]) + Value(multmap_l2_u1_4[all_pair[i]]) + Value(multmap_l2_u2_2[all_pair[i]]) + Value(multmap_l2_u2_3[all_pair[i]]) + Value(multmap_l2_u2_4[all_pair[i]]) + Value(multmap_l2_u3_2[all_pair[i]]) + Value(multmap_l2_u3_3[all_pair[i]]) + Value(multmap_l2_u3_4[all_pair[i]]) + Value(multmap_l2_u4_2[all_pair[i]]) + Value(multmap_l2_u4_3[all_pair[i]]) + Value(multmap_l2_u4_4[all_pair[i]]) + Value(multmap_l2_u5_2[all_pair[i]]) + Value(multmap_l2_u5_3[all_pair[i]]) + Value(multmap_l2_u5_4[all_pair[i]]) + Value(multmap_l3_u1_2[all_pair[i]]) + Value(multmap_l3_u1_3[all_pair[i]]) + Value(multmap_l3_u1_4[all_pair[i]]) + Value(multmap_l3_u2_2[all_pair[i]]) + Value(multmap_l3_u2_3[all_pair[i]]) + Value(multmap_l3_u2_4[all_pair[i]]) + Value(multmap_l3_u3_2[all_pair[i]]) + Value(multmap_l3_u3_3[all_pair[i]]) + Value(multmap_l3_u3_4[all_pair[i]]) + Value(multmap_l3_u4_2[all_pair[i]]) + Value(multmap_l3_u4_3[all_pair[i]]) + Value(multmap_l3_u4_4[all_pair[i]]) + Value(multmap_l3_u5_2[all_pair[i]]) + Value(multmap_l3_u5_3[all_pair[i]]) + Value(multmap_l3_u5_4[all_pair[i]]) + Value(multmap_l4_u1_2[all_pair[i]]) + Value(multmap_l4_u1_3[all_pair[i]]) + Value(multmap_l4_u1_4[all_pair[i]]) + Value(multmap_l4_u2_2[all_pair[i]]) + Value(multmap_l4_u2_3[all_pair[i]]) + Value(multmap_l4_u2_4[all_pair[i]]) + Value(multmap_l4_u3_2[all_pair[i]]) + Value(multmap_l4_u3_3[all_pair[i]]) + Value(multmap_l4_u3_4[all_pair[i]]) + Value(multmap_l4_u4_2[all_pair[i]]) + Value(multmap_l4_u4_3[all_pair[i]]) + Value(multmap_l4_u4_4[all_pair[i]]) + Value(multmap_l4_u5_2[all_pair[i]]) + Value(multmap_l4_u5_3[all_pair[i]]) + Value(multmap_l4_u5_4[all_pair[i]]) + Value(multmap_l5_u1_2[all_pair[i]]) + Value(multmap_l5_u1_3[all_pair[i]]) + Value(multmap_l5_u1_4[all_pair[i]]) + Value(multmap_l5_u2_2[all_pair[i]]) + Value(multmap_l5_u2_3[all_pair[i]]) + Value(multmap_l5_u2_4[all_pair[i]]) + Value(multmap_l5_u3_2[all_pair[i]]) + Value(multmap_l5_u3_3[all_pair[i]]) + Value(multmap_l5_u3_4[all_pair[i]]) + Value(multmap_l5_u4_2[all_pair[i]]) + Value(multmap_l5_u4_3[all_pair[i]]) + Value(multmap_l5_u4_4[all_pair[i]]) + Value(multmap_l5_u5_2[all_pair[i]]) + Value(multmap_l5_u5_3[all_pair[i]]) + Value(multmap_l5_u5_4[all_pair[i]]) + Value(multmap_l6_u1_2[all_pair[i]]) + Value(multmap_l6_u1_3[all_pair[i]]) + Value(multmap_l6_u1_4[all_pair[i]]) + Value(multmap_l6_u2_2[all_pair[i]]) + Value(multmap_l6_u2_3[all_pair[i]]) + Value(multmap_l6_u2_4[all_pair[i]]) + Value(multmap_l6_u3_2[all_pair[i]]) + Value(multmap_l6_u3_3[all_pair[i]]) + Value(multmap_l6_u3_4[all_pair[i]]) + Value(multmap_l6_u4_2[all_pair[i]]) + Value(multmap_l6_u4_3[all_pair[i]]) + Value(multmap_l6_u4_4[all_pair[i]]) + Value(multmap_l6_u5_2[all_pair[i]]) + Value(multmap_l6_u5_3[all_pair[i]]) + Value(multmap_l6_u5_4[all_pair[i]]);
	
		ofs1 << hi + iref + mult << " " << hi << " " << iref << " " << mult << " " ;
		ofs1<<all_pair[i].first<<" "<<all_pair[i].second<<endl;
	} 
	return 0;
}
