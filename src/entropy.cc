#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
#include<iostream>
#include<map>
#include<math.h>
using namespace std;
#define MIN_SIZE  atoi(argv[3])
#define MAX_SIZE  atoi(argv[4])
#define RAND_NUM  atoi(argv[5])
#define RAND atoi(argv[6])
//	../../Stat C Matrix01  2 10  1000000 1 path

void sort(float a[], int n, int b[])
{
	for (int i = 2; i < n + 1; i++)
	{
		for (int j = n; j >= i; j--)
		{
			if (a[j] > a[j - 1])
			{
				float tmp1 = a[j];
				a[j] = a[j - 1];
				a[j - 1] = tmp1;
				int tmp2 = b[j];
				b[j] = b[j - 1];
				b[j - 1] = tmp2;
			}
		}
	}
}

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

float Get_entropy_value(vector<vector<int> > Matrix02);

int main(int argc, char* argv[])
{
	int i = 0, j, k, p, size, size1, size2, size3, size_0, j1, k1, k2, k3, num = 0, num_ILP, v_0, N;
	int SIZE, SIZE1, SIZE2, SIZE3, SIZE4, SIZE5, SIZE6, M_size;
	float q, f, cov;
	vector<int> vec_0, int_line, neighbour_set, pathway, pathway1;
	vector<vector<int> >  Matrix01, Matrix02, C, C_copy, Sub_network;
	vector<vector<float> >  vec, vec1;
	vector<float> line, entropy, entropy1;
	vector<string> Line;
	vector<int> Num, Num_1;
	map<int, int> map1, map2;
	string temp;
	string sto = ",";
	char  szFileName3[100];
	char  szFileName6[100];
	char  szFileName7[100];
	ofstream ofs3, ofs6, ofs7;
	ifstream ifs1;

	ifs1.open(argv[1]);  //influence graph C
	while (getline(ifs1, temp))
	{
		SIZE = temp.size();
		if (SIZE == 0)
		{
			C.push_back(int_line);  //the first line of C is null
		}
		else
		{
			split(temp, sto, Line);
			size = Line.size();
			for (k = 0; k < size; k++)
			{
				int_line.push_back(atoi(Line[k].c_str()));
			}
			C.push_back(int_line);
			int_line.clear();
			Line.clear();
		}
	}
	ifs1.close();
	C_copy = C;
	N = C.size() - 1;
	k = 0;
	for (i = 1; i <= N; i++)
	{
		if (!C_copy[i].empty())
		{
			Sub_network.push_back(int_line);
			Sub_network[k].push_back(i);
			size = Sub_network[k].size();
			Sub_network[k].insert(Sub_network[k].end(), C_copy[i].begin(), C_copy[i].end());
			size1 = Sub_network[k].size();
			while (size < size1)
			{
				for (j = size; j < size1; j++)
				{
					size2 = Sub_network[k].size();
					v_0 = Sub_network[k][j];
					if (!C_copy[v_0].empty())
					{
						size3 = C_copy[v_0].size();
						for (p = 0; p < size3; p++)
						{
							if (C_copy[v_0][p] != i)
							{
								for (q = 0; q < size2; q++)
								{
									if (C_copy[v_0][p] == Sub_network[k][q])
										break;
								}
								if (q == size2)
									Sub_network[k].push_back(C_copy[v_0][p]);
							}

						}
						C_copy[v_0].erase(C_copy[v_0].begin(), C_copy[v_0].end());
					}
				}
				size = size1;
				size1 = Sub_network[k].size();
			}
			k++;  //k: zhe number of subnetwork
		}
	}

	ifstream ifs2;
	ifs2.open(argv[2]);  //Matrix01
	while (getline(ifs2, temp))
	{
		SIZE = temp.size();
		if (SIZE == 0)
		{
			Matrix01.push_back(int_line);
		}
		else
		{
			split(temp, sto, Line);
			size = Line.size();
			for (k = 0; k < size; k++)
			{
				int_line.push_back(atoi(Line[k].c_str()));
			}
			Matrix01.push_back(int_line);
			int_line.clear();
			Line.clear();
		}
	}
	ifs2.close();
	M_size = Matrix01.size();

	srand((unsigned)time(NULL));
	j = C.size() - 1;  // the first line of C is meaningless and j is the gene number
	//get permutation function
	Num.push_back(0);
	for (i = 1; i <= j; i++)
	{
		Num_1.push_back(i);
	}
	for (i = 1; i <= j; i++)
	{
		j1 = Num_1.size();
		k = rand() % j1;
		Num.push_back(Num_1[k]);  //random sequence of genes
		Num_1.erase(Num_1.begin() + k);
	}
	Num_1.clear();

	if (RAND)
	{
		sprintf(szFileName3, "%sEx_entropy.txt", argv[7]);
		ofs3.open(szFileName3);
	}
	if (!RAND)
	{
		sprintf(szFileName6, "%sEx_entropy.txt", argv[7]);
		ofs6.open(szFileName6);
	}
	sprintf(szFileName7, "%sEx_entropy_%d.log", argv[7], RAND_NUM);
	ofs7.open(szFileName7);

	for (SIZE = MIN_SIZE; SIZE <= MAX_SIZE; SIZE++)
	{
		vec_0.clear();
		k1 = 0;   //record number of genes in subnetwork with size larger than SIZE-1
		SIZE1 = Sub_network.size();
		for (i = 0; i < SIZE1; i++)
		{
			SIZE2 = Sub_network[i].size();
			if (SIZE2 >= SIZE)
			{
				vec_0.insert(vec_0.end(), Sub_network[i].begin(), Sub_network[i].end());  //gene, which can be selected
				k1 += SIZE2;
			}
		}
		num_ILP = k1;
		SIZE3 = vec_0.size();
		if (SIZE3 == 0)
		{
			ofs7 << "Error! No network has size greater than " << SIZE << "." << endl;
		}
		else
		{
			num = 0;
			while (num < RAND_NUM)  //select modules 1000000 times randomly
			{
				map1.clear();
				size_0 = vec_0.size();
				k = rand() % size_0;
				v_0 = vec_0[k];
				pathway.push_back(v_0);
				map1[v_0] = 1;
				size = pathway.size();
				neighbour_set.insert(neighbour_set.end(), C[v_0].begin(), C[v_0].end());
				SIZE4 = C[v_0].size();
				for (i = 0; i < SIZE4; i++)
				{
					map1[C[v_0][i]] = 1;
				}
				while (size < SIZE)  //get modules(pathway) with size of SIZE
				{
					size_0 = neighbour_set.size();
					k = rand() % size_0;
					v_0 = neighbour_set[k];
					pathway.push_back(v_0);
					size++;
					if (size == SIZE)
						break;
					neighbour_set.erase(neighbour_set.begin() + k);
					SIZE5 = C[v_0].size();
					for (i = 0; i < SIZE5; i++)
					{
						if (!map1[C[v_0][i]])
						{
							neighbour_set.push_back(C[v_0][i]);
							map1[C[v_0][i]] = 1;
						}
					}
				}
				num++;
				neighbour_set.clear();
				if (RAND)
				{
					cov = 0;
					for (j = 0; j < M_size; j++)  //for each patient
					{
						k3 = 0;
						for (i = 0; i < SIZE; i++)
						{
							k = pathway[i];
							k2 = Matrix01[j][k];
							int_line.push_back(k2);
							k3 += k2;
						}
						if (k3 > 0)
						{
							cov = cov + 1;
							Matrix02.push_back(int_line);
						}
						int_line.clear();
					}
					f = Get_entropy_value(Matrix02);
					entropy.push_back(f);
					Matrix02.clear();
				}
				if (!RAND)
				{
					for (i = 0; i < SIZE; i++)
					{
						j = pathway[i];
						pathway1.push_back(Num[j]);
					}
					cov = 0;
					for (j = 0; j < M_size; j++)
					{
						k3 = 0;
						SIZE6 = pathway1.size();
						for (i = 0; i < SIZE6; i++)
						{
							k = pathway1[i];
							k2 = Matrix01[j][k];
							int_line.push_back(k2);
							k3 += k2;
						}
						if (k3 > 0)
						{
							cov++;
							Matrix02.push_back(int_line);
						}
						int_line.clear();
					}
					f = Get_entropy_value(Matrix02);
					entropy1.push_back(f);
					Matrix02.clear();
					pathway1.clear();
				}
				pathway.clear();
			}
		}
		if (RAND)
		{
			sort(entropy.begin(), entropy.end(), greater<float>());

			ofs3 << SIZE;
			ofs7 << SIZE;
			size = entropy.size();
			for (k = 0; k < size; k++)
			{
				ofs3 << " " << entropy[k];
				if (k == (int)(size * 0.05) || k == (int)(size * 0.1) || k == (int)(size * 0.15) || k == (int)(size * 0.2) || k == (int)(size * 0.25) || k == (int)(size * 0.3))
					ofs7 << " " << entropy[k];
			}
			ofs3 << endl;
			ofs7 << endl;
			entropy.clear();
		}
		if (!RAND)
		{
			sort(entropy1.begin(), entropy1.end(), greater<float>());
			ofs6 << SIZE;
			ofs7 << SIZE;
			size = entropy1.size();
			for (k = 0; k < size; k++)
			{
				ofs6 << " " << entropy1[k];
				if (k == (int)(size * 0.05) || k == (int)(size * 0.1) || k == (int)(size * 0.15) || k == (int)(size * 0.2) || k == (int)(size * 0.25) || k == (int)(size * 0.3))
					ofs7 << " " << entropy1[k];
			}
			ofs6 << endl;
			ofs7 << endl;
			entropy1.clear();
		}
	}
	return 0;
}


float Get_entropy_value(vector<vector<int> > Matrix02)
{
	int i, j, row, col, num, num1;
	double p;
	float p1;

	row = Matrix02.size();
	col = Matrix02[0].size();
	int* uniq = new int[col + 1];
	for (j = 0; j < col; j++)
		uniq[j] = 0;
	num1 = 0;
	for (i = 0; i < row; i++)
	{
		num = 0;
		for (j = 0; j < col; j++)
		{
			num += Matrix02[i][j];
		}
		if(num==1)
			for (j = 0; j < col; j++)
			{
				if (Matrix02[i][j] == 1)
					uniq[j]++;
			}
		if (num)
			num1++;
	}
	p = p1 = 0;
	for (j = 0; j < col; j++)
	{
		p = (double)uniq[j] / num1;
		if(p != 0)
			p1 -= (p * log(p) / log(col));
	}
	delete []uniq;
	return p1;
}