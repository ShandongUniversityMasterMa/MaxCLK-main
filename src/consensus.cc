#include<fstream>
#include<algorithm>
#include<vector>
#include<string>
#include<iostream>
#include<stdlib.h>
#include<map>
using namespace std;
#define Weight_B atoi(argv[4])
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

int main (int argc, char *argv[])
{
	vector<string> vec, vec1, vec2, vec3;
	vector<string> Line, Line1;
	string temp, temp1;
	Line.clear();
	Line1.clear();
	ifstream ifs, ifs1, ifs2;
	ofstream ofs1, ofs2;
	string sto = " ";
	vector<int> int_line;
	map<string,int> mapboth;
	map<string,int> maponly;
	map<string, int> map1, map2, map3, map4, map6, map7, map10, mapd;
	map<int, int> map11;
	map<pair<int, string>, int> map8, map9;
	map<pair<string, string>, int> map5;
	vector<pair<string, string> > pair2, pair3, all_pair, pair1;
	vector<pair<int, string> > pair1_1, pair1_2, pair2_1, pair2_2, pairlinker;
	vector<vector<string> >  Matrix, CORE, CORE1, CORE2, CORE3, CORE4, CORE5, VEC;
	vector<vector<int> > Linker;
	vector<string> linker, W3, W2, W2n, W1, W1n, W13;
    
	ifs.open(argv[5]);  
	while (getline(ifs, temp))
	{
		mapd[temp] = 1;  
	}
	ifs.close();
	ifs.open(argv[1]);
	ofs1.open(argv[2]);
	ofs2.open(argv[3]);
	int i = 0, j = 0, j1, j1n, j2, j2n, j3, j4, k = 0, k1, i1, i2, i3, i4, i5, i6, kk, index, indictor, SIZE, SIZE1, SIZE2;

	while (getline(ifs, temp))  
	{
		split(temp, sto, Line);
		i1 = atoi(Line[1].c_str());
		i2 = atoi(Line[2].c_str());
		i3 = atoi(Line[3].c_str());
		if (i1 > 0) { i4 = 1; }
		else { i4 = 0; }
		if (i2 > 0) { i5 = 1; }
		else { i5 = 0; }
		if (i3 > 0) { i6 = 1; }
		else { i6 = 0; }
		if (i4 + i5 + i6 == 3)
		{
			all_pair.push_back(make_pair(Line[4], Line[5]));
			if (!map2[Line[4]])
			{
				W3.push_back(Line[4]);  
				map2[Line[4]] = 1;  
			}
			if (!map2[Line[5]])
			{
				W3.push_back(Line[5]);
				map2[Line[5]] = 1;
			}
		}
		Line.clear();
	}
	ifs.close();

	/*Combine genes in all_pair and get CORE network*/
	kk = 0;
	index = 0;
	indictor = 0;
	SIZE = all_pair.size();  
	while (kk < SIZE)  
	{
		if (indictor == 0)
		{
			index++;
			map1[all_pair[kk].first] = index;
			map1[all_pair[kk].second] = index;  
			map5[all_pair[kk]] = index;  
		}
		indictor = 0;
		for (i = kk; i < SIZE; i++)
		{
			if (map1[all_pair[i].first] > 0 && map1[all_pair[i].second] > 0)
			{
				map5[all_pair[i]] = map1[all_pair[i].first];
			}

			if (map1[all_pair[i].first] > 0 && map1[all_pair[i].second] == 0)
			{
				map1[all_pair[i].second] = map1[all_pair[i].first];
				map5[all_pair[i]] = map1[all_pair[i].first];
				indictor = 1;
			}
			if (map1[all_pair[i].first] == 0 && map1[all_pair[i].second] > 0)
			{
				map1[all_pair[i].first] = map1[all_pair[i].second];
				map5[all_pair[i]] = map1[all_pair[i].first];
				indictor = 1;
			}
		}
		for (i = 0; i < SIZE; i++)
		{
			if (map5[all_pair[i]] == 0)
			{
				break;
			}
		}
		kk = i;
	}

	for (i = 1; i <= index; i++)  
	{
		CORE.push_back(Line);  
		for (j = 0; j < SIZE; j++)
		{
			if (map5[all_pair[j]] == i)
			{
				if (!map7[all_pair[j].first])
				{
					CORE[i - 1].push_back(all_pair[j].first);
					map7[all_pair[j].first] = i;
				}
				if (!map7[all_pair[j].second])
				{
					CORE[i - 1].push_back(all_pair[j].second);
					map7[all_pair[j].second] = i;
				}
			}
		}
	}
	map5.clear();
	map7.clear();
	all_pair.clear();

	SIZE = CORE.size();
	if (SIZE > 0)
	{
		for (i = 0; i < SIZE; i++)
		{
			SIZE1 = CORE[i].size();
			if (SIZE1 > 0)
			{
				sort(CORE[i].begin(), CORE[i].end());
			}
		}
		for (i = 0; i < SIZE; i++)
		{
			CORE[i].push_back("+");  
		}
	}

	CORE2 = CORE;  
	ifs.open(argv[1]);
	while (getline(ifs, temp))  
	{
		split(temp, sto, Line);
		i1 = atoi(Line[1].c_str());
		i2 = atoi(Line[2].c_str());
		i3 = atoi(Line[3].c_str());
		if (i1 > 0) { i4 = 1; }
		else { i4 = 0; }
		if (i2 > 0) { i5 = 1; }
		else { i5 = 0; }
		if (i3 > 0) { i6 = 1; }
		else { i6 = 0; }

		if (i4 + i5 + i6 == 2 && map2[Line[4]] + map2[Line[5]] != 2)
		{
			if (map2[Line[4]])  
			{
				j = map1[Line[4]] - 1;
				if (!map8[make_pair(j, Line[5])])
				{
					pair2_1.push_back(make_pair(j, Line[5]));  
					if (!map3[Line[5]])
					{
						W2.push_back(Line[5]);  
						map3[Line[5]] = 1;  
					}
					map8[make_pair(j, Line[5])] = 1;
				}
			}
			if (map2[Line[5]])
			{
				j = map1[Line[5]] - 1;
				if (!map8[make_pair(j, Line[4])])
				{
					pair2_1.push_back(make_pair(j, Line[4]));
					if (!map3[Line[4]])
					{
						W2.push_back(Line[4]);
						map3[Line[4]] = 1;
					}
					map8[make_pair(j, Line[4])] = 1;
				}
			}
			if (!map2[Line[4]] && !map2[Line[5]])  
			{
				pair3.push_back(make_pair(Line[4], Line[5]));  
			}
			/* Line[4],Line[5] have been sequenced */
		}
		Line.clear();
	}
	ifs.close();

	SIZE = pair2_1.size();
	for (i = 0; i < SIZE; i++)
	{
		map4[pair2_1[i].second]++;
	}
	for (i = 0; i < SIZE; i++)
	{
		if (map4[pair2_1[i].second] == 1)  
		{
			j = pair2_1[i].first;
			CORE[j].push_back(pair2_1[i].second);
			map7[pair2_1[i].second] = j + 1;  
		}
	}
	SIZE = CORE.size();
	if (SIZE > 0)
	{
		for (i = 0; i < SIZE; i++)
		{
			SIZE1 = CORE[i].size();
			for (j = 0; j < SIZE1; j++)
			{
				if (CORE[i][j] == "+")
				{
					break;
				}
			}
			if (j != SIZE1 - 1)
			{
				sort(CORE[i].begin() + j + 1, CORE[i].end());
			}
			CORE[i].push_back("*");  
		}
	}

	CORE3 = CORE;
	SIZE = pair2_1.size();
	for (i = 0; i < SIZE; i++)
	{
		if (map4[pair2_1[i].second] > 1)  
		{
			j = pair2_1[i].first;
			CORE[j].push_back(pair2_1[i].second);
			if (!map6[pair2_1[i].second])
			{
				linker.push_back(pair2_1[i].second);
				Linker.push_back(int_line);  
				SIZE1 = Linker.size();
				map6[pair2_1[i].second] = SIZE1;
				k = SIZE1 - 1;
				Linker[k].push_back(j + 1);
			}
			else
			{
				k = map6[pair2_1[i].second] - 1;
				Linker[k].push_back(j + 1);
			}
		}
	}
	map4.clear();
	SIZE = Linker.size();
	for (i = 0; i < SIZE; i++)
	{
		sort(Linker[i].begin(), Linker[i].end());
	}
	SIZE = CORE.size();
	if (SIZE > 0)
	{
		for (i = 0; i < SIZE; i++)
		{
			SIZE1 = CORE[i].size();
			for (j = 0; j < SIZE1; j++)
			{
				if (CORE[i][j] == "*")
				{
					break;
				}
			}
			if (j != SIZE1 - 1)
			{
				sort(CORE[i].begin() + j + 1, CORE[i].end());
			}
			CORE[i].push_back("-");  
		}
	}
	CORE4 = CORE;

	ifs.open(argv[1]);
	while (getline(ifs, temp))
	{
		split(temp, sto, Line);
		i1 = atoi(Line[1].c_str());
		i2 = atoi(Line[2].c_str());
		i3 = atoi(Line[3].c_str());
		if (i1 > 0) { i4 = 1; }
		else { i4 = 0; }
		if (i2 > 0) { i5 = 1; }
		else { i5 = 0; }
		if (i3 > 0) { i6 = 1; }
		else { i6 = 0; }
		
		if (i4 + i5 + i6 == 1 && map2[Line[4]] + map2[Line[5]] == 1)  
		{
			if (map2[Line[4]] && !map3[Line[5]])  
			{
				j = map1[Line[4]] - 1;
				if (!map8[make_pair(j, Line[5])])
				{
					pair1_1.push_back(make_pair(j, Line[5]));  
					map8[make_pair(j, Line[5])] = 1;
				}
			}
			if (map2[Line[5]] && !map3[Line[4]])
			{
				j = map1[Line[5]] - 1;
				if (!map8[make_pair(j, Line[4])])
				{
					pair1_1.push_back(make_pair(j, Line[4]));
					map8[make_pair(j, Line[4])] = 1;

				}
			}
		}
		Line.clear();
	}
	ifs.close();
	map4.clear();
	SIZE = pair1_1.size();
	for (i = 0; i < SIZE; i++)
	{
		map4[pair1_1[i].second]++;
	}
	for (i = 0; i < SIZE; i++)
	{
		if (map4[pair1_1[i].second] == 1)
		{
			pair1_2.push_back(pair1_1[i]);  
			j = pair1_1[i].first;
			CORE[j].push_back(pair1_1[i].second);
			W1.push_back(pair1_1[i].second);  
		}
	}
	map4.clear();
	SIZE = CORE.size();
	if (SIZE > 0)
	{
		for (i = 0; i < SIZE; i++)
		{
			SIZE1 = CORE[i].size();
			for (j = 0; j < SIZE1; j++)
			{
				if (CORE[i][j] == "-")
				{
					break;
				}
			}
			if (j != SIZE1 - 1)
			{
				sort(CORE[i].begin() + j + 1, CORE[i].end());
			}
			CORE[i].push_back("--");
		}
	}
	CORE5 = CORE;
	
	kk = 0;
	index = 0;
	indictor = 0;
	SIZE = pair3.size();
	while (kk < SIZE)  
	{
		if (indictor == 0)
		{
			index++;
			map4[pair3[kk].first] = index;
			map4[pair3[kk].second] = index;  
			map5[pair3[kk]] = index;  
		}
		indictor = 0;
		for (i = kk; i < SIZE; i++)
		{
			if (map4[pair3[i].first] > 0 && map4[pair3[i].second] > 0)
			{
				map5[pair3[i]] = map4[pair3[i].first];
			}

			if (map4[pair3[i].first] > 0 && map4[pair3[i].second] == 0)
			{
				map4[pair3[i].second] = map4[pair3[i].first];
				map5[pair3[i]] = map4[pair3[i].first];
				indictor = 1;
			}
			if (map4[pair3[i].first] == 0 && map4[pair3[i].second] > 0)
			{
				map4[pair3[i].first] = map4[pair3[i].second];
				map5[pair3[i]] = map4[pair3[i].first];
				indictor = 1;
			}
		}
		for (i = 0; i < SIZE; i++)
		{
			if (map5[pair3[i]] == 0)
			{
				break;
			}
		}
		kk = i;

	}
	map4.clear();

	for (i = 1; i <= index; i++)
	{
		CORE1.push_back(Line);
		for (j = 0; j < SIZE; j++)
		{
			if (map5[pair3[j]] == i)
			{
				if (!map4[pair3[j].first])
				{
					CORE1[i - 1].push_back(pair3[j].first);
					map4[pair3[j].first] = i;
				}
				if (!map4[pair3[j].second])
				{
					CORE1[i - 1].push_back(pair3[j].second);
					map4[pair3[j].second] = i;
				}
			}
		}
	}
	map4.clear();
	
    SIZE=CORE.size();
	for (i = 0; i < SIZE; i++)
	{
		SIZE1 = CORE[i].size();
		for (j = 0; j < SIZE1; j++)
		{
			ofs1 << CORE[i][j] << " ";
		}
		ofs1 << endl;
	}
    ofs1<<endl;

	SIZE = CORE1.size();  //CORE1: type-2
	if (SIZE > 0)
	{
		for (i = 0; i < SIZE; i++)
		{
			sort(CORE1[i].begin(), CORE1[i].end());
			kk = 0;
			SIZE1 = CORE1[i].size();
			for (j = 0; j < SIZE1; j++)
			{
				if (map3[CORE1[i][j]])
				{
					kk++;  
				}
			}
			if (kk > 0)
			{
				for (j = 0; j < SIZE1; j++)
				{
					ofs1 << CORE1[i][j];
					if (map3[CORE1[i][j]])  
					{
						if (map7[CORE1[i][j]])  
						{
							ofs1 << "(" << map7[CORE1[i][j]] << ") " << " ";
						}
						else  
						{
							k = map6[CORE1[i][j]] - 1;
							ofs1 << "(";
							SIZE2 = Linker[k].size() - 1;
							for (k1 = 0; k1 < SIZE2; k1++)
							{
								ofs1 << Linker[k][k1] << ",";
							}
							ofs1 << Linker[k][k1] << ") " << " ";
						}
					}
					else
					{
						ofs1 << " ";
					}
				}
				ofs1 << endl;
			}
			else
			{
				for (j = 0; j < SIZE1; j++)
				{
					ofs1 << CORE1[i][j] << " ";
				}
				ofs1 << endl;
			}
		}
	}
	map3.clear();
	map4.clear();

	j3 = 0;
	SIZE = W3.size();  
	for (i = 0; i < SIZE; i++)
	{
		map3[W3[i]] = 1;  
		if (mapd[W3[i]])  
		{
			j3++;  
		}
	}
	cout << "Method A: type-1 weight 3: " << j3 << "/" << SIZE << "=" << (double)j3 / SIZE << " ";

	j2 = 0;
	SIZE = W2.size();
	for (i = 0; i < SIZE; i++)
	{
		map3[W2[i]] = 1;  
		if (mapd[W2[i]])
		{
			j2++;  
		}
	}
	cout << "type-1 weight 2: " << j2 << "/" << SIZE << "=" << (double)j2 / SIZE << " ";

	SIZE = pair3.size();  
	for (i = 0; i < SIZE; i++)
	{
		if (!map3[pair3[i].first])
		{
			W2n.push_back(pair3[i].first);  
			map3[pair3[i].first] = 1;
		}
		if (!map3[pair3[i].second])
		{
			W2n.push_back(pair3[i].second);
			map3[pair3[i].second] = 1;  
		}
	}
	j2n = 0;
	SIZE = W2n.size();
	for (i = 0; i < SIZE; i++)
	{
		map3[W2n[i]] = 1;
		if (mapd[W2n[i]])
		{
			j2n++;  
		}
	}
	cout << "type-2 weight2 not in type-1: " << j2n << "/" << SIZE << "=" << (double)j2n / SIZE << " ";

	j4 = 0;
	SIZE = W1.size();  
	for (i = 0; i < SIZE; i++)
	{
		if (mapd[W1[i]])
		{
			j4++;  
		}
	}
	j1 = 0;
	k1 = 0;
	for (i = 0; i < SIZE; i++)
	{
		if (!map3[W1[i]])
		{
			map3[W1[i]] = 1;  
			if (mapd[W1[i]])
			{
				j1++;  
			}
			k1++;
		}
	}

	ifs.open(argv[1]);
	while (getline(ifs, temp))
	{
		split(temp, sto, Line);
		if (!map3[Line[4]])
		{
			W1n.push_back(Line[4]);  
			map3[Line[4]] = 1;
		}
		if (!map3[Line[5]])
		{
			W1n.push_back(Line[5]);
			map3[Line[5]] = 1;
		}
		Line.clear();
	}
	ifs.close();

	j1n = 0;
	SIZE = W1n.size();
	for (i = 0; i < SIZE; i++)
	{
		if (mapd[W1n[i]])
		{
			j1n++;  
		}
	}
	SIZE = W3.size();
	SIZE1 = W2.size();
	SIZE2 = W2n.size();
	cout << "type-1 type-2: " << j3 + j2 + j2n << "/" << SIZE + SIZE1 + SIZE2 << "=" << (double)(j3 + j2 + j2n) / (SIZE + SIZE1 + SIZE2) << " ";


	i = 0;
	map2.clear();
	map3.clear();
	vec.clear();
	ifs.open(argv[1]);
	while (getline(ifs, temp))
	{
		split(temp, sto, Line);
		if (atoi(Line[0].c_str()) >= Weight_B)
		{
			if (!map2[Line[4]] && !map2[Line[5]])
			{
				i++;
				VEC.push_back(Line1);
				vec.push_back(Line[4]);
				vec.push_back(Line[5]);
				VEC[i - 1].push_back(Line[4]);
				VEC[i - 1].push_back(Line[5]);
				map2[Line[4]] = i;
				map2[Line[5]] = i;
			}
			else if (!map2[Line[4]] && map2[Line[5]])
			{
				k = map2[Line[5]];
				VEC[k - 1].push_back(Line[4]);
				vec.push_back(Line[4]);
				map2[Line[4]] = k;
			}
			else if (map2[Line[4]] && !map2[Line[5]])
			{
				k = map2[Line[4]];
				VEC[k - 1].push_back(Line[5]);
				vec.push_back(Line[5]);
				map2[Line[5]] = k;
			}
			else
			{
				if (map2[Line[4]] != map2[Line[5]])
				{
					if (map2[Line[4]] < map2[Line[5]])
					{
						j = map2[Line[4]];
						k = map2[Line[5]];
					}
					else
					{
						j = map2[Line[5]];
						k = map2[Line[4]];
					}
					map11[k] = 1;
					SIZE = VEC[k - 1].size();
					for (j1 = 0; j1 < SIZE; j1++)
					{
						VEC[j - 1].push_back(VEC[k - 1][j1]);
						map2[VEC[k - 1][j1]] = j;
					}
				}
			}
		}
		Line.clear();
	}
	ifs.close();
	SIZE = VEC.size();
	for (i = 0; i < SIZE; i++)
	{
		if (!map11[i + 1])
		{
			sort(VEC[i].begin(), VEC[i].end());
			SIZE1 = VEC[i].size();
			for (k = 0; k < SIZE1; k++)
			{
				ofs2 << VEC[i][k] << " ";
			}
			ofs2 << endl;
		}
	}
	sort(vec.begin(), vec.end());
	k = 0;
	SIZE = vec.size();
	for (i = 0; i < SIZE; i++)
	{
		if (mapd[vec[i]])
		{
			k++;
		}
	}
	cout << "Method B: " << k << "/" << SIZE << "=" << (double)k / SIZE << endl;

    return 0;
}