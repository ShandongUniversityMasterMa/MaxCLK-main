#include<fstream>
#include<vector>
#include<string>
#include<iostream>
#include<stdlib.h>
#include<algorithm>
#include<cmath>
#include<eigen3/Eigen/Dense>
#include<math.h>
using namespace Eigen;
using namespace std;

#define nSet atoi(argv[2])
#define min_mutation_num atoi(argv[3])
#define errEx atof(argv[4])
#define errCov atof(argv[5])

const int maxn = 350;
int nV; //nV vertex number, nEdge edge number
vector<int> edge[maxn];  
vector<int> Vertex;  
vector<vector<int>> tmcliques;  
MatrixXi GenGramG;

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

//int Partition(vector<int>& arr, int indx[], int start, int end);
void swap(vector<int>& arr, int i, int j);
void swap2(int indx[], int i, int j);
void QuickSort(vector<int>& arr, int indx[], int left, int right);

void dfs(int pos, int nset)  
{

    int j, SIZE, SIZE1;
    if(Vertex.size() == nset)
    {
       tmcliques.push_back(Vertex);  
       return;
    }
    
    SIZE = edge[pos].size();
    for(int i = 0; i < SIZE; i++){
        SIZE1 = Vertex.size();
        for (j = 0; j < SIZE1; j++){
            if (!GenGramG(edge[pos][i], Vertex[j]))
                break;
        }  
        if(j==Vertex.size()){
            Vertex.push_back(edge[pos][i]);
            dfs(edge[pos][i],nset);
            Vertex.pop_back();  
        }
    }
}

// ./maxCliques mutMat.txt nSet min_mutation_num  errEx errCov  Out_Dir/module.txt
//             vector<vector<int> >& mat, int nSet, int min_mutation_num, int errEx, float errCov,  Out_Dir/module.txt

int main(int argc, char* argv[])
{
    int k = 0, col = 0, ncol, nCol, row, nrow, nRow, nMat, nsize, nvsum, ntmp;
    vector<int> int_line,colSum, Genes, tGenes, tCov, Cov, uniq;
    vector<float> tentropy, entropy, Exc, tExc, CovEx, tCovEx;
    vector<string> Line;
    vector<vector<int> > tmutMat, clkMat, mcliques,module;
    int temp1, temp2;
    float f, f1, f2, f3;
	string temp;
    string sto=",";
    ifstream ifs;
	ofstream ofs;


    ifs.open(argv[1]); 
	while(getline(ifs, temp))
	{
		k++;
		if(k==1)
		{
			split(temp, sto, Line);
			nMat=atoi(Line[col].c_str());  
		}
		else
		if(k==2)
		{
			split(temp, sto, Line);
			nsize=Line.size();
			for(col=0; col<nsize; col++)
			{
              tGenes.push_back(atoi(Line[col].c_str()));  
			}
		}
		else
		{
			split(temp, sto, Line);
			nsize=Line.size();
			for(col=0; col<nsize; col++)
			{
              int_line.push_back(atoi(Line[col].c_str()));
			}
	        tmutMat.push_back(int_line);  
			int_line.clear();
		}
		Line.clear();
	}
	ifs.close();
    nRow = tmutMat.size();  
    nCol = tmutMat[0].size();

    //sum column  and then sort the column of matrix according to the sort colSum

	int* indx = new int[nCol];  
    nvsum = 0; 
    for (col = 0; col < nCol; col++){
        temp1 = 0;
        for (row = 0; row < nRow; row++){
            temp1 = temp1 + tmutMat[row][col]; 
        }
        colSum.push_back(temp1);  
        if(colSum[col] >= min_mutation_num)
            nvsum++;  
    }
	// BubbleSort(colSum, nCol, indx, false); //false : sort the indx but first one.
    for (int m = 0; m < nCol; m++)
		indx[m] = m;	    
    QuickSort(colSum, indx, 1, nCol - 1); //sort the indx but first one.
    for (col = 0; col < nCol; col++)
        Genes.push_back(tGenes[indx[col]]);  

    MatrixXi mutMat = MatrixXi::Zero(1,1);
    if(nvsum == col)
    { 
        mutMat.resize(nRow, nvsum);    //include the first one although it's sum maybe less min_mutation_num
        for (row = 0; row < nRow; row++)
        {
            for (col = 0; col < nvsum; col++)
            {
                mutMat(row,col) = tmutMat[row][indx[col]];
            }
        }
    }
    else{
        nvsum++;
        mutMat.resize(nRow, nvsum);    //include the first one although it's sum maybe less min_mutation_num
        for (row = 0; row < nRow; row++)
        {
            for (col = 0; col < nvsum; col++)
            {
                mutMat(row,col) = tmutMat[row][indx[col]];
            }
        }
    }
   
    // Calculate genes Gramian matrix
    MatrixXi GenGram = mutMat.transpose() * mutMat;
    nV = GenGram.rows();  

    // build Genes Gramian Graph
    for (int i = 0 ; i < nV ; i++)
    {
        edge[i].clear();  
    }
    GenGramG.setZero(nV,nV);  
    for(nrow = 0; nrow < nV; nrow++)
    {
        for(ncol = nrow + 1; ncol < nV; ncol++)
        {
            if(GenGram(nrow, ncol) <= errEx)
            {
               GenGramG(nrow,ncol) = 1;
               GenGramG(ncol,nrow) = 1;
               edge[nrow].push_back(ncol);  
            }
        }
    }

//-------------------------------------------------------------
    // Run maximal cliques algorithm.
    tmcliques.clear();
    Vertex.push_back(0);
    dfs(0, nSet);  

    tCov.clear();
    tentropy.clear();
    nsize = tmcliques.size();
    for(row =0; row < nsize; row++)
    {
        temp1 = 0;
        for (col = 0; col < nSet; col++)
            uniq.push_back(0);
        for (int i = 0; i < nRow; i++)
        {
            temp2 = 0;
            for (col = 0; col < nSet; col++)
            {
                k = tmcliques[row][col];
                temp2 += mutMat(i, k);
            }
            if (temp2)
                temp1++;
            if(temp2 == 1)
                for (col = 0; col < nSet; col++)
                {
                    k = tmcliques[row][col];
                    if (mutMat(i, k) == 1)
                    {
                        uniq[col]++;  
                        break;
                    }
                }
        }
        f = 0;
        f1 = 0;
        f2 = 0;
        for (col = 0; col < nSet; col++)
        {
            k = tmcliques[row][col];
            f = (float)uniq[col] / temp1;
            if(f)
                f1 -= (float)(f * log(f) / log(nSet));
            f2 += ((float)uniq[col] / colSum[k]);
        }
        f2 = f2 / nSet;
        tCov.push_back(temp1);  
        tentropy.push_back(f1);  
        tExc.push_back(f2);  //Ex
        ntmp = tCov[row];
        f3 = ntmp * f2 / nMat;
        tCovEx.push_back(f3);  //CovEx
        uniq.clear();
    }
    
//-------------------------------------------------------------------
    // get the original label for each maximal clique
    mcliques.clear();
    int_line.clear();
    for(row = 0; row < nsize; row++)
    {
        for(col = 0; col < nSet ; col++ )
        {
            k = tmcliques[row][col];
            ntmp=tGenes[indx[k]];
            int_line.push_back(ntmp);
        }
        mcliques.push_back(int_line);
        int_line.clear();
    }
    tmcliques.clear();
    Vertex.clear();
       
    for(row = 0; row < nsize; row++)  
    {
        if(tCov[row] > errCov * (float)nMat)  
        {
            module.push_back(mcliques[row]);
            Cov.push_back(tCov[row]);
            entropy.push_back(tentropy[row]);
            Exc.push_back(tExc[row]);
            CovEx.push_back(tCovEx[row]);
        }
    }
    nsize = module.size();
    mcliques.clear();
    tCov.clear();
    tentropy.clear();
    tExc.clear();
    tCovEx.clear();

//--------------------------------------------------------------
// CLKmodule_nSet.txt:  s Coverage Cov Ex CovEx tHnat genelist, to create SKmodule_nSet.txt (= Smodule_nSet.txt) by filter_CLK.cc

    ofs.open(argv[6],ios_base::app);  
    for(row = 0; row < nsize; row++)
    {
        ofs << module[row][0] << " " << Cov[row] << " " << (double)Cov[row] / nMat << " " << Exc[row] << " " << CovEx[row] << " " << entropy[row] << " ";
        for(col = 0; col < nSet; col++ )
        {
            ofs << module[row][col] << " ";
        }
        ofs<<endl;
    }
    ofs.close();
    return 0;
}


void swap(vector<int>& arr, int i, int j)
{
    int temp;
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}
void swap2(int indx[], int i, int j)
{
    int temp;
    temp = indx[i];
    indx[i] = indx[j];
    indx[j] = temp;
}

void QuickSort(vector<int>& arr, int indx[], int left, int right)
{
    int i, last;
    void swap(vector<int>&arr, int i, int j);
    void swap2(int indx[], int i, int j);
    if (left >= right)
        return;
    swap(arr, left, (left + right) / 2);
    swap2(indx, left, (left + right) / 2);
    last = left;
    for (i = left + 1; i <= right; i++)
    {
        if (arr[i] > arr[left])
        {
            swap(arr, ++last, i);
            swap2(indx, last, i);
        }
    }
    swap(arr, left, last);
    swap2(indx, left, last);
    QuickSort(arr, indx, left, last - 1);
    QuickSort(arr, indx, last + 1, right);
}
