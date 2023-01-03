#include<ctime>
#include<fstream>
#include<vector>
#include<string>
#include <algorithm>
#include<iostream>
#include<map>
using namespace std;

//../minset_greedy  Matrix01  Result.txt(SKmodule_2.txt)  Considered_genes CSKmodule_2.txt

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
	int i=0,j=0,i1,j1,k,p,size,pat,M_size,M_size0,new_M_size,kk=0;
	int sum=0,sum1=0,sum2=0,indictor,SIZE,SIZE1,SIZE2;
	float q;

    vector<int> Num,int_line,Indictor,Patients1,Patients2,Patients3,Sum;
	vector<string> Line;
	vector<float> line;
	vector<vector<float> > Group;
	vector<vector<int> > Matrix01,Matrix02,Matrix_gene,Gene_group;

	map<int,int> map1;   //To distinguish selected groups from non-selected groups
	map<int,string> map2;
	map<string,int> map5;
	map<int,int> map3;   //To distinguish selected genes from non-selected genes
	map<int,int> map4;   //To get genes in gene groups(results)

	ofstream ofs,ofs1;
	ofs.open(argv[4]);
	ofs1.open(argv[5]);

        string temp;
	string sto=",";
	string sto1=" ";

	ifstream ifs;
	ifs.open(argv[3]);
	while(getline(ifs, temp))
	{
		split(temp, sto1,  Line) ;
		i=atoi(Line[1].c_str());
		map2[i]=Line[0];
		map5[Line[0]]=i;
		Line.clear();
	}
	ifs.close();

	ifs.open(argv[2]);
	while(getline(ifs, temp))
	{
		split(temp, sto1,  Line) ;
		size=Line.size();
		for(k=6;k<size;k++)
		{
			i=map5[Line[k]];
			int_line.push_back(i);
			map4[i]=1;
		}
		Gene_group.push_back(int_line);
		int_line.clear();
		for(k=0;k<6;k++)
		{
			q=atof(Line[k].c_str());
			line.push_back(q);
		}
		Group.push_back(line);
		line.clear();
		Line.clear();
	}
	ifs.close();

        ifs.open(argv[1]);
	while(getline(ifs, temp))
	{
		if(temp!="")
		{
			split(temp, sto,  Line) ;
			size=Line.size();
			for(k=0;k<size;k++)
			{
				i=atoi(Line[k].c_str());
				int_line.push_back(i);
			}
	         	Matrix01.push_back(int_line);
			int_line.clear();
			Line.clear();
		}
	}
	ifs.close();
	M_size=Matrix01.size();
	M_size0=Matrix01[0].size();
	//ofs<<"Patient number covered by the considered genes which are also in the PPI network: "<<M_size<<endl;

	for(i=1;i<M_size0;i++)
	{
		if(!map4[i])
		for(j=0;j<M_size;j++)
		{
			Matrix01[j][i]=0;
		}
	}
	for(i=0;i<M_size;i++)
	{
		p=0;
		for(j=0;j<M_size0;j++)
		{
			p+=Matrix01[i][j];
		}
		if(p>0)
		{
			Matrix02.push_back(Matrix01[i]);    // patients having mutations in gene groups are considered only
		}
	}
	Matrix01.clear();

	new_M_size=Matrix02.size();
	for(i=0;i<new_M_size;i++)
	{
		Patients3.push_back(0);
	}

	Matrix_gene.push_back(int_line);
	for(i=1;i<M_size0;i++)
	{
		for(j=0;j<new_M_size;j++)
		{
			if(Matrix02[j][i]==1)
			{
				int_line.push_back(j);
			}
		}
		Matrix_gene.push_back(int_line);
		int_line.clear();
	}
	//ofs<<"Patient number covered by the significant modules: "<<new_M_size<<endl;


	sum=0;
	sum1=0;
	kk=0;
	size=Gene_group.size();
	while(sum<new_M_size)
	{
		kk++;
		for(i=0;i<size;i++)
		{
			if(!map1[i])
			{
				Patients1=Patients3;
				SIZE=Gene_group[i].size();
				for(j=0;j<SIZE;j++)
				{
					i1=Gene_group[i][j];
					if(!map3[i1])
					{
						SIZE1=Matrix_gene[i1].size();
						for(j1=0;j1<SIZE1;j1++)
						{
							pat=Matrix_gene[i1][j1];
							Patients1[pat]=1;
						}
					}
				}
				sum2=0;
				SIZE2=Patients1.size();
				for(j=0;j<SIZE2;j++)
				{
					sum2+=Patients1[j];
				}
				if(sum2>sum1)
				{
					sum1=sum2;
					indictor=i;
					Patients2=Patients1;
				}
			}
		}
		Indictor.push_back(indictor);
		SIZE=Gene_group[indictor].size();
		for(i=0;i<SIZE;i++)
		{
			j=Gene_group[indictor][i];
			map3[j]=1;
		}
		int_line.push_back(sum1-sum);
		sum=sum1;
		Patients3=Patients2;
		map1[indictor]=1;
		Sum.push_back(sum);
	}

	SIZE=Indictor.size();
	for(i=0;i<SIZE;i++)
	{
		k=Indictor[i];
		//ofs<<int_line[i]<<" "<<Sum[i]<<" "<<(double)Sum[i]/new_M_size;
		ofs<<int_line[i];
		SIZE1=Group[k].size();
		for(j=0;j<SIZE1;j++)
		{
			ofs<<" "<<Group[k][j];
		}
		SIZE2=Gene_group[k].size();
		for(j=0;j<SIZE2;j++)
		{
			ofs<<" "<<map2[Gene_group[k][j]];
			if(j==0)
			{
				ofs1<<map2[Gene_group[k][j]];
			}
			else
			{
				ofs1<<" "<<map2[Gene_group[k][j]];
			}
		}
		ofs<<endl;
		ofs1<<endl;
	}
	ofs.close();
	ofs1.close();


	return 0;
}
