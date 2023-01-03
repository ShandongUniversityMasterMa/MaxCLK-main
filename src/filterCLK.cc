#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
#include<iostream>
#include<map>
using namespace std;
//	../../filter_CLK CLKmodule_i.txt Considered_genes  rand_stat_covexc.txt 0.8 0.05
#define Value_Ex atof(argv[4])
#define P_Exentropy atof(argv[5])

void sort(float a[], int n, int b[])
{
	for (int i = 1; i < n + 1; i++)
		for (int j = n; j >= i; j--)
			if (a[j] < a[j - 1])
			{
				float tmp1 = a[j];
				a[j] = a[j - 1];
				a[j - 1] = tmp1;
				int tmp2 = b[j];
				b[j] = b[j - 1];
				b[j - 1] = tmp2;
			}
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


int main (int argc, char *argv[])
{
	int i = 0, j, k, nsize, size2, n;
	float f, f1, f2, q;
	vector<int> int_line1, int_line2;
	vector<vector<int>> GENE, VALUE;
	vector<vector<float>> STAT_Ex_entropy, Group;
    vector<float> line, Exc, Ex_entropy;
	vector<string> Line;
	map<int,int> map1;
	map<int,string> map2;
    string temp;
	string sto1=" ";
	ifstream ifs;

    ifs.open(argv[1]);
	while(getline(ifs, temp))
	{
		split(temp, sto1,  Line) ;
		nsize=Line.size();
		for(k=0;k<2;k++)
		{
			int_line1.push_back(atoi(Line[k].c_str()));
		}
		for(k=1;k<6;k++)
		{
			q = atof(Line[k].c_str());
			line.push_back(q);
		}
		for(k=6;k<nsize-1;k++)   //the last character is space
		{
			int_line2.push_back(atoi(Line[k].c_str()));
		}
		sort(int_line2.begin(), int_line2.end());
	    GENE.push_back(int_line2);   // genelist
		VALUE.push_back(int_line1);  // s Coverage
		Group.push_back(line);       //   Coverage Cov Ex CovEx tHnat
		int_line1.clear();
		int_line2.clear();
		line.clear();
		Line.clear();
	}
	ifs.close();
	ifs.open(argv[2]);  
	while(getline(ifs, temp))
	{
		split(temp, sto1, Line) ;
		map2[atoi(Line[1].c_str())] = Line[0];
		Line.clear();
	}
	ifs.close();

	ifs.open(argv[3]);
	i=0;
	while(getline(ifs, temp))
	{
		split(temp, sto1, Line);
		nsize = Line.size();
		for(k=1;k<nsize;k++)
		{
			line.push_back(atof(Line[k].c_str()));
		}
	    STAT_Ex_entropy.push_back(line);
		j=atoi(Line[0].c_str());
		map1[j]=i;
		i++;
		line.clear();
		Line.clear();
	}
	ifs.close();

	size2 = GENE.size();
	n = size2 - 1;
	float* a = new float[size2 + 1];
	int* b = new int[size2 + 1];
	for (i = 0; i <= n; i++)
	{
		b[i] = i;
	}
	
	for(i=0;i<size2;i++)
	{
		nsize=GENE[i].size(); //   GENElist
        // Group:  Coverage Cov Ex CovEx Exentropy
		f = Group[i][2];
		Exc.push_back(f);
		f1 = Group[i][4];
		Ex_entropy.push_back(f1);	
		k = map1[nsize];
		nsize = STAT_Ex_entropy[k].size();
		for (j = 0; j < nsize; j++)
		{
			if (f1 + 0.000001 > STAT_Ex_entropy[k][j])
			{
				break;
			}
		}
		f2 = (double)j / nsize;
		a[i]=f2;
//		std::cout << "f1==="<<f1<< "f2="<< f2 << "STAT_COVEXC[k][0]" <<STAT_COVEXC[k][0]<< endl; //for test
	}
	sort(a, n, b);
	ofstream ofs;
	ofs.open(argv[6]);

// SKmodule_nSet.txt:  Coverage Cov Ex CovEx entropy p_value_entropy  genelist
	for (i = 0; i <= n; i++)
	{
		k = b[i];
		if (Exc[k] >= Value_Ex && a[i] <= P_Exentropy)
		{
			for (j = 0; j < 5; j++)
			{
				ofs << Group[k][j] << " ";
			}
			ofs << a[i] << " ";
			nsize = GENE[k].size();
			for (j = 0; j < nsize - 1; j++)
			{
				ofs << map2[GENE[k][j]] << " ";
			}
			ofs << map2[GENE[k][j]] << endl;
		}
		if (a[i] > P_Exentropy)
			break;
	}
    return 0;
}
