#include<fstream>
#include<vector>
#include<string>
#include<iostream>
#include<stdlib.h>
#include<algorithm>
#include<map>
using namespace std;

#define D atof(argv[4])

//./prep_data mutation_data.txt  PPI_matrix.txt  PPI_index_file.txt   15(delta)  Out_Dir genes_considered
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
	int nlist,edge_num,ch,ch1=0;
	int i=0,j=0,k=0,nMat,nGenes;
	int nsize=0,SIZE;
    int row, col,nRow, nCol;
	float DELTA;
	vector<vector<int>> C, Matrix01, snvs, snvs1;
	vector<vector<int>> mutated_genes;
	vector<int> indictor, Gene_network;
	vector<int> int_line, Genes, isolate_genes;
	vector<float> InfluenceMatrix, Inf;
	vector<string> Line;
	map<string, int> map1, map2;
	map<int,int> map3;
	map<int, string> map5;
	string temp, sto = ",", sto1 = " ";
	int_line.clear();
	Line.clear();
	char szFileName[100];
	char szFileName3[100];
	ifstream ifs, ifs3;
	ofstream ofs1, ofs3;

    ifs.open(argv[6]); 
	ch = ifs.get();
	if (ch != EOF)
	{
		ch1 = 1;
		ifs3.open(argv[6]);
		while (getline(ifs3, temp))
		{
			map1[temp] = 1;  
		}
		ifs3.close();
	}
    ifs.close();
	ifs.open(argv[3]);  
	i = 0;
	while (getline(ifs, temp))
	{
		i++;
		split(temp, sto1, Line);
		k = atoi(Line[0].c_str());
		map2[Line[1]] = k;  
		map5[k] = Line[1];  
		Line.clear();
	}
	nlist = i;  
	ifs.close();

	ifs.open(argv[1]);  
	if (ch1 == 1)
	{
		while (getline(ifs, temp))
		{
			split(temp, sto1, Line);
			SIZE = Line.size();
			for (i = 1; i < SIZE; i++)  
			{
				k = map2[Line[i]];
				if (map1[Line[i]] && k)
				{
					int_line.push_back(k);
				}
			}
			snvs.push_back(int_line);  
			Line.clear();
			int_line.clear();
		}
	}
	else
	{
		while (getline(ifs, temp))
		{
			split(temp, sto1, Line);
			nsize = Line.size();
			for (i = 1; i < nsize; i++)
			{
				k = map2[Line[i]];
				if (k)
				{
					int_line.push_back(k);
				}
			}
			snvs.push_back(int_line);
			Line.clear();
			int_line.clear();
		}
	}
	ifs.close();
	map1.clear();

	nRow = snvs.size();  
	for (row = 0; row < nRow; row++)
	{
		if (!snvs[row].empty())
		{
			snvs1.push_back(snvs[row]); 
			nCol = snvs[row].size();
			for (col = 0; col < nCol; col++)
			{
				if (!map3[snvs[row][col]])
				{
					Genes.push_back(snvs[row][col]); 
					map3[snvs[row][col]] = 1;  
				}
			}
		}
	}
	map3.clear();

	sort(Genes.begin(),Genes.end());
	nsize=Genes.size();  
	for (i = 0; i < nsize; i++)
	{
		k = Genes[i];
		map3[k] = i + 1;  
	}
	nMat = snvs1.size();  
	nGenes = Genes.size(); 
	sprintf(szFileName3, "%sWrite.log", argv[5]);  
	ofs3.open(szFileName3);
	ofs3 << "Number of all patients: " << snvs.size() << endl;
	ofs3 << "Number of covered patients " << nMat << endl;
	ofs3 << "Number of considered genes: " << nGenes << endl;
	sprintf(szFileName, "%sGene_num.txt", argv[5]);
	ofs1.open(szFileName);
	ofs1 << nGenes;
	ofs1.close();

	for (col = 0; col <= nGenes; col++)
	{
		int_line.push_back(0);    
	}
	for (row = 0; row < nMat; row++)
	{
		Matrix01.push_back(int_line);
		nsize = snvs1[row].size();
		for (k = 0; k < nsize; k++)
		{
			j = snvs1[row][k];
			col = map3[j];  
			Matrix01[row][col] = 1; 
		}
	}
	int_line.clear();
	ifs.close();

	ifs.open(argv[2]); 
	while (getline(ifs, temp))
	{
		split(temp, sto, Line);
		nsize = Line.size();  
		for (j = 0; j < nsize; j++)
		{
			InfluenceMatrix.push_back(atof(Line[j].c_str()));  
		}
		Line.clear();
	}
	ifs.close();
	// cout<< "size of PPI matrix: "<< nsize << endl;

	for (i = 0; i < nGenes - 1; i++)
	{
		for (j = i + 1; j < nGenes; j++)
		{
			row = Genes[i] - 1;  
			col = Genes[j] - 1;  
			nRow = row * nlist + col;
			nCol = col * nlist + row;
			if (InfluenceMatrix[nRow] >= InfluenceMatrix[nCol])  
			{
				Inf.push_back(InfluenceMatrix[nCol]);
			}
			else
			{
				Inf.push_back(InfluenceMatrix[nRow]);
			}
		}
	}
	sort(Inf.rbegin(), Inf.rend());
	j = int(D * nGenes / 2 - 1);
	DELTA = Inf[j];  

	ofs3 << "delta: " << DELTA << endl;
	for (i = 0; i <= nGenes; i++)
	{
		C.push_back(int_line);  
	}
	edge_num = 0;
	for (i = 0; i < nGenes - 1; i++)
	{
		for (j = i + 1; j < nGenes; j++)
		{
			row = Genes[i] - 1;
			col = Genes[j] - 1;
			nRow = row * nlist + col;
			nCol = col * nlist + row;
			if (InfluenceMatrix[nRow] >= DELTA && InfluenceMatrix[nCol] >= DELTA)  
			{
				C[i + 1].push_back(j + 1);
				C[j + 1].push_back(i + 1);
				edge_num++;  
			}
		}
	}
	InfluenceMatrix.clear();
	ofs3 << "edge number: " << edge_num << endl;
	ofs3.close();

//--------------------------------------------------------------------------------
	sprintf(szFileName, "%sC", argv[5]);
	ofs1.open(szFileName);
	nRow = C.size();
	for (i = 0; i < nRow; i++)  
	{
		nCol = C[i].size();
		if (nCol == 0)
			ofs1 << endl;
		else
		{
			if (nCol == 1)
				ofs1 << C[i][0] << endl;
			else
			{
				for (j = 0; j < nCol - 1; j++)
				{
					ofs1 << C[i][j] << ",";
				}
				ofs1 << C[i][j] << endl;
			}
		}
	}
	ofs1.close();

	sprintf(szFileName, "%sMatrix01", argv[5]);
	ofs1.open(szFileName);
	for (i = 0; i < nMat; i++)
	{
		nCol = Matrix01[i].size();
		if (nCol == 0)
			ofs1 << endl;
		else
		{
			if (nCol == 1)
				ofs1 << Matrix01[i][0] << endl;
			else
			{
				for (j = 0; j < nCol - 1; j++)
				{
					ofs1 << Matrix01[i][j] << ",";
				}
				ofs1 << Matrix01[i][j] << endl;
			}
		}
	}
	ofs1.close();

	sprintf(szFileName, "%sConsidered_genes", argv[5]);
	ofs1.open(szFileName);
	nsize = Genes.size();
	for (i = 0; i < nsize; i++)
	{
		k = Genes[i];
		ofs1 << map5[k] << " " << i + 1 << endl;  
	}
	ofs1.close();

	return 0;
}
