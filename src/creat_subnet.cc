#include<fstream>
#include<vector>
#include<string>
#include<iostream>
#include<stdlib.h>
#include<algorithm>
#include<map>
using namespace std;

#define NUM atof(argv[3])
#define MIN_T atoi(argv[4])
#define MAX_T atoi(argv[5])

inline static void split(std::string src, std::string token, vector<std::string>& vect)
{
	int nend = 0;
	int nbegin = 0;
	while (nend != -1)
	{
		nend = src.find_first_of(token, nbegin);
		if (nend == -1)
			vect.push_back(src.substr(nbegin, src.length() - nbegin));
		else
			vect.push_back(src.substr(nbegin, nend - nbegin));
		nbegin = nend + 1;
	}
}

//./creat_mMat  C  Matrix01  300(local_network_size)  2(min_size) 5(max_size) Out_Dir

int main (int argc, char *argv[])
{
	int nMAX = 2;
	int i = 0, j = 0, k = 0;
	int nlist, nMat, nGenes, nInd, p, q;
	int nsize = 0, nset, nGnet, nnGnet, s, SIZE;
	int row, nRow, col, nCol;
	vector<vector<int> > C, Matrix01, subMat;
	vector<vector<int> > mutated_genes;
	vector< int > indictor, Gene_network;
	vector< int > int_line, Genes, isolate_genes;
	vector< float > InfluenceMatrix, Inf;
	vector<string> Line;

//---------------------------------------------------------//
	map<int, int> map3;
	string temp, sto = ",", sto1 = " ";
	int_line.clear();
	char  szFileName1[100];
	char  szFileName3[100];
	ifstream ifs, ifs3;
	ofstream ofs1, ofs3;

//------------------------------------------------------------------//
    ifs.open(argv[1]);  
	while (getline(ifs, temp))
	{
		if (temp.size() == 0)
		{
			C.push_back(int_line);
		}
		else
		{
			split(temp, sto, Line);
			nsize = Line.size();
			for (k = 0; k < nsize; k++)
			{
				int_line.push_back(atoi(Line[k].c_str()));
			}
			C.push_back(int_line);
			int_line.clear();
			Line.clear();
		}
	}
	ifs.close();
	nGenes = C.size();

    ifs.open(argv[2]);  
	while (getline(ifs, temp))
	{
		if (temp.size() == 0)
		{
			Matrix01.push_back(int_line);
		}
		else
		{
			split(temp, sto, Line);
			nsize = Line.size();
			for (k = 0; k < nsize; k++)
			{
				int_line.push_back(atoi(Line[k].c_str()));
			}
			Matrix01.push_back(int_line);
			int_line.clear();
			Line.clear();
		}
	}
	ifs.close();
	nMat = Matrix01.size();

	sprintf(szFileName3,"%sWrite.log",argv[6]);
	ofs3.open(szFileName3,ios::app);
	k = 0;
	nlist = nGenes;  
	//cout << nlist << endl;
	//cout << C[1][0] << endl;
  	// cout << "s==" <<N<< endl;
//------------------------------------------------------------------//
	for (s = 1; s < nlist; s++)
	{
		if (C[s].empty())  
		{
			isolate_genes.push_back(s);
		}
		else
		{
			Gene_network.clear();
			indictor.clear();
			Gene_network.push_back(s);
			Gene_network.insert(Gene_network.end(), C[s].begin(), C[s].end());
			ofs3 << "Search radius information of gene " << s << ":" << endl;

			indictor.push_back(1);
			SIZE = C[s].size();
			if (SIZE + 1 > NUM)
			{
				indictor.push_back(SIZE + 1);
				ofs3 << "1," << SIZE + 1 << endl;
			}
			else
			{
				nsize = 1;
				nGnet = Gene_network.size();
				while (nsize < nGnet && nGnet <= NUM)
				{
					indictor.push_back(nGnet);
					for (k = nsize; k < nGnet; k++)
					{
						nnGnet = Gene_network.size();
						p = Gene_network[k];
						nInd = C[p].size();
						for (q = 0; q < nInd; q++)
						{
							for (i = 0; i < nnGnet; i++)
							{
								if (C[p][q] == Gene_network[i])
								{
									break;
								}
							}
							if (i == nnGnet)
							{
								Gene_network.push_back(C[p][q]);
							}
						}
					}
					nsize = nGnet;
					nGnet = Gene_network.size();
				}
				nInd = indictor.size();
				for (i = 0; i < nInd; i++)
				{
					ofs3 << indictor[i] << ",";
				}
				ofs3 << nGnet << endl;
			}

		    //--------------------------------------------------------------//
			if (nsize >= MAX_T)
			{
				nMAX = MAX_T;
			}
			else if (nsize >= MIN_T)
			{
				nMAX = nsize;
			}
			if (nsize >= MIN_T)
			{
				for (nset = MIN_T; nset <= nMAX; nset++)
				{
					SIZE = indictor.size();
					if (SIZE >= nset)
					{
						nInd = indictor[nset - 1] - 1;
					}
					else
					{
						nInd = indictor.back() - 1;
					}
					Genes.clear();
					for (k = 0; k <= nInd; k++)  
					{
						nGenes = Gene_network[k];
						Genes.push_back(nGenes);
					}

					map3.clear();
					sort(Genes.begin() + 1, Genes.end());  
					nGnet = Genes.size();
					for (i = 0; i < nGnet; i++)
					{
						k = Genes[i];
						map3[k] = i;  
					}

					mutated_genes.clear();
					for (j = 0; j < nMat; j++)
					{
						for (k = 0; k <= nInd; k++)
						{
							nGenes = Gene_network[k];
							if (Matrix01[j][nGenes] == 1)
							{
								int_line.push_back(nGenes);
							}
						}
						if (!int_line.empty())
						{
							mutated_genes.push_back(int_line);
							int_line.clear();
						}
					}
					int_line.clear();
					nRow = mutated_genes.size();  

					vector<vector<int>> subMat(nRow, vector<int>(nInd + 1, 0));  

					for (row = 0; row < nRow; row++)
					{
						nsize = mutated_genes[row].size();
						for (col = 0; col < nsize; col++)
						{
							nCol = mutated_genes[row][col];
							i = map3[nCol];
							subMat[row][i] = 1;  
						}
					}
					int_line.clear();

					//-------------------------------------------------------------------
					sprintf(szFileName1, "%smMat%d_t%d.txt", argv[6], s, nset);
					ofs1.open(szFileName1);
					ofs1 << nMat << "," << endl;
					ofs1 << s;
					for (j = 1; j <= nInd; j++)
					{
						ofs1 << "," << Genes[j];
					}
					ofs1 << endl;
					for (i = 0; i < nRow; i++)
					{
						nCol = subMat[i].size();
						if (nCol == 0)
							ofs1 << endl;
						else
						{
							if (nCol == 1)
								ofs1 << subMat[i][0] << endl;
							else
							{
								for (j = 0; j < nCol - 1; j++)
								{
									ofs1 << subMat[i][j] << ",";
								}
								ofs1 << subMat[i][j] << endl;
							}
						}
					}
					ofs1.close();
				}
			}
		}
	}
	ofs3.close();
	return 0;
}
