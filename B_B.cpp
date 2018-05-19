// Naive algorithm for building suffix array of a given text
#include <iostream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <time.h>
#include <ctime>
#include <sstream>
#include <unordered_map>
using namespace std;

// Structure to store information of a suffix
struct suffix
{
	int index;
	char *suff;
};
struct factor
{
  char *str;
  vector<int> position;

};
 struct node
{
string pattern;
int mul;

int distance;

};
 struct node2
{
string pattern;
int mul;
int num;
int distance;

};
string LCS="";
string LCS_dummy="";
string LCS1="";
string check;
vector<int> indices;
vector<int> LCS_index;
vector<string> LCS_result;
 ofstream myfile;
int index_L=0;
int bestsolution=0;
 int v=0,temp_v=0;
vector<string> stringg;
 unordered_map<string, int> umap;
string f_;
vector<int> marker1,marker2,marker3;
vector<int> oc_num1,oc_num2;
 char *txt[100];
  int j=0;
  vector<node> NodeID;
int count_DFS=0;
vector< vector<int> > vec;
vector<string > result;
 clock_t t1,t2;
 int dis=0,dis1=0;
 string key="";
// thread th1;
// A comparison function used by sort() to compare two suffixes
int cmp(struct suffix a, struct suffix b)
{
	return strcmp(a.suff, b.suff) < 0? 1 : 0;
}
void  DFS( vector< vector<node> > Node,node Node_element,string lcs,int dis,int dis1,vector <int> c);
void DFS_part(vector< vector<node> > Node,string LCS,string LCS_dummy,int bestsolution,int a,int dis,int count_item,vector <int> marker1,vector <int> marker_dummy);
// This is the main function that takes a string 'txt' of size n as an
// argument, builds and return the suffix array for the given string
int *buildSuffixArray(char *txt, int n)
{
	// A structure to store suffixes and their indexes
	struct suffix suffixes[n];

	// Store suffixes and their indexes in an array of structures.
	// The structure is needed to sort the suffixes alphabatically
	// and maintain their old indexes while sorting
	for (int i = 0; i < n; i++)
	{
		suffixes[i].index = i;
		suffixes[i].suff = (txt+i);


	}

	// Sort the suffixes using the comparison function
	// defined above.
	sort(suffixes, suffixes+n, cmp);

	// Store indexes of all sorted suffixes in the suffix array
	int *suffixArr = new int[n];
	for (int i = 0; i < n; i++)
	{
		suffixArr[i] = suffixes[i].index;

      //cout<<suffixes[i].index<<'\n';
		//cout<<suffixes[i].suff<<'\n';
	}


	// Return the suffix array
	return suffixArr;
}

int find_occurrance_num(char *txt, int n,char *pat,int occ_num)
{
	// A structure to store suffixes and their indexes
	struct suffix suffixes[n];

	// Store suffixes and their indexes in an array of structures.
	// The structure is needed to sort the suffixes alphabatically
	// and maintain their old indexes while sorting
	for (int i = 0; i < n; i++)
	{
string s=&txt[0];
char *str = new char[1];
s= s.substr(i,1);
strcpy(str,s.c_str());
		suffixes[i].index = i;
		suffixes[i].suff =str;

//cout<< str<<'\n';
	}

	// Sort the suffixes using the comparison function
	// defined above.
	sort(suffixes, suffixes+n, cmp);

	// Store indexes of all sorted suffixes in the suffix array
	int *suffixArr = new int[n];
	for (int i = 0; i < n; i++)
	{
		suffixArr[i] = suffixes[i].index;
//cout<<suffixes[i].index<<'\n';
		//cout<<suffixes[i].suff<<'\n';

	}

int m = strlen(pat);
int r=0; int i=0;
int index=0,index_low=0;
int count_pat=0;
int mid=0,mid_low=0;
vector<int> intermediate;
	if((n%2)==0)
	r=n;
	else
	r=n-1;
	int l = 0;
	while (l <= r)
	{

		mid = l + (r-l)/2;
		int res = strncmp(pat,txt+suffixArr[mid],m);

		if (res == 0)
		{

      index=suffixArr[mid];
      index_low=suffixArr[mid];
      mid_low=mid-1;

      break;

		}


		if (res < 0) r = mid - 1;

		else l =mid +1;
	}
int *suffixArr2 = new int[2];

while(strncmp(pat,suffixes[mid].suff,1)==0)
{
intermediate.push_back(suffixes[mid].index);
mid++;

}

while(strncmp(pat,suffixes[mid_low].suff,1)==0)
{
intermediate.push_back(suffixes[mid_low].index);
mid_low--;


}

sort(intermediate.begin(),intermediate.end());

return intermediate[occ_num];
}


int* find_occurrance(char *txt, int n,char *pat)
{
	// A structure to store suffixes and their indexes
	struct suffix suffixes[n];

	// Store suffixes and their indexes in an array of structures.
	// The structure is needed to sort the suffixes alphabatically
	// and maintain their old indexes while sorting
	for (int i = 0; i < n; i++)
	{
string s=&txt[0];
char *str = new char[1];
s= s.substr(i,1);
strcpy(str,s.c_str());
		suffixes[i].index = i;
		suffixes[i].suff =str;

//cout<< str<<'\n';
	}

	// Sort the suffixes using the comparison function
	// defined above.
	sort(suffixes, suffixes+n, cmp);

	// Store indexes of all sorted suffixes in the suffix array
	int *suffixArr = new int[n];
	for (int i = 0; i < n; i++)
	{
		suffixArr[i] = suffixes[i].index;
//cout<<suffixes[i].index<<'\n';
		//cout<<suffixes[i].suff<<'\n';

	}

int m = strlen(pat);
int r=0; int i=0;
int index=0,index_low=0;
int count_pat=0;
int mid=0,mid_low=0;

	if((n%2)==0)
	r=n;
	else
	r=n-1;
	int l = 0;
	while (l <= r)
	{

		mid = l + (r-l)/2;
		int res = strncmp(pat,txt+suffixArr[mid],m);

		if (res == 0)
		{

      index=suffixArr[mid];
      index_low=suffixArr[mid];
      mid_low=mid-1;

      break;

		}


		if (res < 0) r = mid - 1;

		else l =mid +1;
	}
int *suffixArr2 = new int[2];

while(strncmp(pat,suffixes[mid].suff,1)==0)
{
if(suffixes[mid].index<=index)
{
index=suffixes[mid].index;


}
count_pat++;
mid++;

}

while(strncmp(pat,suffixes[mid_low].suff,1)==0)
{
if(suffixes[mid_low].index<index_low)
{
 //cout<<suffixes[mid_low].index<<'\n';
index_low=suffixes[mid_low].index;


}
count_pat++;
mid_low--;


}

if(index_low<=index)
{
index=index_low;

}
else
{
index=index;
}
suffixArr2[0]=n-1-(index+1);

 suffixArr2[1]=count_pat;
if(count_pat==0)
{
suffixArr2[0]=0;

 suffixArr2[1]=0;
}
 count_pat=0;
index=0;
return suffixArr2;
}

 node make_Node(const string &s,const int a, const int b)
    {
    node nodes;
    nodes.pattern = s;
    nodes.mul = a;

    nodes.distance = b;
        return nodes;
    }
set<char> str_intersection(string s[],int k)
{
char *str1;
set<char> s1;
set<char> s2;
set<char> s3;
for(int i=0;i<s[0].size();i++)
{

s1.insert(s[0][i]);



}
for(int i=0;i<s[1].size();i++)
{

s2.insert(s[1][i]);



}

  vector<char> v_intersection;
set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), back_inserter(v_intersection));
s2.clear();

s2.insert(v_intersection.begin(),v_intersection.end());
v_intersection.clear();

    s1.clear();
    for(int l=2;l<k;l++)
    {

for(int i=0;i<s[l].size();i++)
{

s1.insert(s[l][i]);

}

  set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), back_inserter(v_intersection));
  s2.clear();
  s2.insert(v_intersection.begin(),v_intersection.end());
  v_intersection.clear();
  s1.clear();
    }

return s2;


}

void show(  vector< vector<node> > Node)
{


for(int i=0; i<Node.size();i++)
{
    cout<<'\n';
    for(int l=0;l<Node[i].size();l++)

   {

   cout<<Node[i][l].pattern<<" "<<Node[i][l].mul<<" "<<Node[i][l].distance<<" ";
   cout<<'\t';
    }
}

}
vector<int> check_position_validity(string s,string t,string s1,vector<int> marker){
 marker2.clear();string sk;
 vector<int> m=marker;
 int k=0,l=0;

for(int kl=0;kl<j;kl++)
{
string copy=&txt[kl][0];
int n=(copy.substr(m[kl],copy.size()).find(s));
int m_=(copy.substr(m[kl],copy.size()).find(t));
sk=copy.substr(m[kl],copy.size());
//cout<<sk<<'\n';
//if(s=="X"){cout<<"markX:"<<copy.substr(marker1[kl],copy.size())<<'\n';}
if(n==-1)
{
//cout<<s<<" "<<t<<" "<<sk;
k=-1;
break;

// return 0;
}
if(n==-1 && m_>=0 )
{
k=-1;
break;
//cout<<s<<" "<<t<<" "<<sk;
 //return 0;
}
if((n==-1)||(sk.size()==0))
{
k=-1;

//cout<<sk<<'\n';
break;
 //return 0;


}
//else
{

{ marker2.push_back(sk.find(s)+1);

if(sk.size()==0)
{
k=-1;
break;

}
//int c=marker1[kl];

  m[kl]+=marker2[kl];
if(s=="X"&&t=="C")
{
 //cout<<sk<<m<<marker2[kl]<<s1<<'\n';

}
l=1;
  }




}

}
 if(k==-1)
 {
m.clear();
 return m;
 }

 else{

  return m;
 }
}
int find_element(vector< vector<node> > Node, string A,int s,int t)
{
int count_flag=0;
//cout<<s<<'\n';
//cout<<t<<'\n';

for(int i=0;i<Node.size();i++)
{
if((s==0))
{

if((Node[i][0].pattern==A)&&(Node[i][0].mul==1))
{
count_flag=i;


break;
}
//else
{
//continue;
}
}
 if((Node[i][0].pattern==A)&&(Node[i][0].distance<=t))
{



count_flag=i;


break;
}
count_flag=i;
if((Node[i][0].pattern==A)&&(Node[i][0].mul==1)&&(Node[i][0].distance>t))
{

count_flag=i;
//continue;
break;

}
//if((Node[i][0].pattern!=A))
{

//count_flag=i;
//continue;
//break;

}
}
//cout<<"index"<<count_flag<<'\n';
return count_flag;
}

int max(int a,int b)
{
//cout<<a<<'\t'<<b<<'\n';
if(a>b)
return a;
else
return b;
}
void DFS_part(vector< vector<node> > Node,string LCS,string LCS_dummy,int bestsolution,int a,int dis,int count_item,vector <int> marker1,vector <int> marker_dummy)

{






}
void  DFS( vector< vector<node> > Node,node Node_element,string lcs,int dis,int dis1,vector <int> c)
{


//cout<<Node_element.pattern<<'\t'<<Node_element.mul<<'\t'<<Node_element.distance<<'\n';
string LCS=lcs;
//cout<<lcs;

string LCS_dummy=lcs;

int dis2=dis;
int dis3;
vector<int> marker1=c,marker_dummy;
marker_dummy.clear();
for(int kl=0;kl<c.size();kl++)
{
//marker1.push_back(c[kl]);
//cout<<marker1[kl]<<'\n';
}
int count_char=0,count_item=0;

//dis=Node_element.distance;
//if(dis>0)
count_item=find_element(Node,Node_element.pattern,dis2,dis1);
//cout<<count_item;
temp_v=0;
LCS_dummy=LCS;

//cout <<"dmy: "+LCS_dummy<<'\n';
  dis3=Node[count_item][0].distance;
//cout<<dis3<<'\n';

//cout<<LCS_dummy<<'\n';
for(int kl=0;kl<c.size();kl++)
{
marker_dummy.push_back(marker1[kl]);
}




//cout<<"found"<<'\n';



//cout<<LCS_dummy.length()<<'\t'<<umap[key]<<'\t'<<bestsolution<<'\n';
for(int j=1;j<Node[count_item].size();j++)
{
v=0;

//cout<<LCS_dummy<<'\t'<<Node[count_item][j].pattern<<'\n';

marker1=check_position_validity(Node[count_item][j].pattern,LCS_dummy.substr(LCS_dummy.size()-1,LCS_dummy.size()),LCS_dummy,marker_dummy);
if(marker1.size()>0)
{
stringstream ss;
ss<<Node[count_item][j].mul;
key=Node[count_item][j].pattern+ss.str();
if (umap.find(key) != umap.end())
{
//cout<<"old"<<'\n';

//cout<<key<<'\t'<<umap[key]<<'\t'<<bestsolution<<'\n';

if(LCS1.length()+umap[key]>bestsolution)
{


LCS=LCS_dummy+Node[count_item][j].pattern;

DFS(Node,Node[count_item][j],LCS,Node[count_item][j].distance,dis3,marker1);



}


}
else if (umap.find(key) == umap.end())
{

//cout<<"new"<<'\n';
LCS=LCS_dummy+Node[count_item][j].pattern;

DFS(Node,Node[count_item][j],LCS,Node[count_item][j].distance,dis3,marker1);

}



}



}

//else
{
//cout<<LCS_dummy.length()<<'\t'<<umap[key]<<'\t'<<bestsolution<<'\n';


}
//cout<<LCS<<'\n';
stringstream ss;
ss<<Node_element.mul;
key=Node_element.pattern+ss.str();
if (umap.find(key) != umap.end())
{
if(LCS.length()>bestsolution)
        {

       // cout<<bestsolution<<'\n';
        bestsolution=LCS.length();
        LCS1=LCS;

        t2=clock();
cout<<LCS1<<'\t'<<LCS1.size()<<'\n';


float diff ((float)t2-(float)t1);
cout<<diff / (double) CLOCKS_PER_SEC<<endl;

        }

}


//cout<<Node_element.pattern<<Node_element.mul<<'\n';


if (umap.find(key) == umap.end())
{
if(LCS1.length()==LCS_dummy.length())
{
umap[key]=1;

}
else
{

umap[key]=LCS1.length()-LCS_dummy.length()+1;

}
//cout<<key<<": "<<umap[key]<<'\n';

}
else if(LCS.length()-LCS_dummy.length()>umap[key])
{


umap[key]=LCS.length()-LCS_dummy.length()+1;
//cout<<key<<": "<<umap[key]<<'\n';



}
else
{


//cout<<key<<": "<<umap[key]<<'\n';
}


}




vector<string> getUniques(vector<string> collection)
{
    vector<string> uniques;
    for (vector<string>::iterator it=collection.begin(); it!=collection.end(); ++it)
    {
        if (find(uniques.begin(), uniques.end(), *it) == uniques.end())
            uniques.push_back(*it);
    }

    return uniques;
}
int main()
{

   t1=clock();





   ifstream file("/home/user/Desktop/motif_discovery/m_1");
    myfile.open ("/home/user/Desktop/motif_discovery/m_2");
   string str;
   string genome_con;
   string genome_concat;
   char* genome_concat_char;
   string substrm[1000],substrn[1000];
   int  k=0,m=0,n=0;
   while ( getline(file, str))
    {

       txt[j] = new char[str.length() + 1];
       strcpy(txt[j], str.c_str());


vec.push_back(vector<int>());

cout<<txt[j]<<"\n";
string genome(txt[j]);
genome_con=genome;
genome_concat= genome_concat+txt[j];
        j++;

    }
for (int i = 0; i < j; i++) {
    vector<int> row; // Create an empty row

    vec.push_back(row); // Add the row to the main vector
}
 vector<int*> txt2;
    genome_concat_char=&genome_concat[0];
    int combined_length=genome_concat.size();
    //cout<<'\n'<<genome_concat_char<<combined_length;

    char *pat[genome_con.size()+1] ;


        int s = genome_con.size();

        int i=0,Z=0;

 vector<factor> sub[j+1];
 vector<string> patterns;
 set<string> pattern;
 set<char *> final_pattern;



vector<string> patterns_new;
//set<string> pattern_new;
set<string> pattern_new_set1;
set<string> pattern_new_set2;
string s1;
string carry="N";
int Mul_sum=0;int fla=0;
for (int i = 0; i < strlen(txt[0]); i++)
{
 string s1=&txt[0][0];
s1= s1.substr(i,1);


 patterns_new.push_back(s1);

}
  vector<string> pattern_new = getUniques(patterns_new);
//patterns_new.erase( unique( patterns_new.begin(), patterns_new.end() ), patterns_new.end() );
   //pattern_new.insert(patterns_new.begin(), patterns_new.end());
cout<<pattern_new.size();
 	int dis=100000,mul=100000;
 for (vector<string>::iterator it=pattern_new.begin(); it!=pattern_new.end(); ++it)
 {
 for(int i=0;i<j;i++)
 {
int ntt=strlen(txt[i]);
char *str = new char[1];
string s2=*it;
//cout<<s2<<'\n';
strcpy(str,s2.c_str());
int * arr=find_occurrance(txt[i], ntt+1,str);
if(arr[1]==0)
{
carry=*it;
}
//cout<<arr[0]<<'\t'<<arr[1]<<'\n';
if((arr[1]<=mul))
{
dis=arr[0];


mul=arr[1];

}

}
//cout<<dis<<'\t'<<mul<<'\n';
Mul_sum=Mul_sum+mul;
if(*it!=carry)
NodeID.push_back(make_Node(*it,mul,dis));
cout<<" ";
    dis=10000000;mul=1000000000;

}

     vector< vector<node> > Node(Mul_sum);
     vector<node> Node_element;
     Node.push_back(Node_element);

vector<node> Nodes;
patterns_new.clear();
 dis=10000000;mul=10000000;
string store_str[j],store_str1[j];
 int * arr;
 int count=0;
 int multiplicity=0;
 int best_length=100000;
for(int i=0;i<NodeID.size();i++)
{
//
char *str = new char[1];
 strcpy(str,NodeID[i].pattern.c_str());
for (int l=0;l<NodeID[i].mul;l++)
{
multiplicity=NodeID[i].mul;
//cout<<NodeID[i].mul<<'\n';
//cout<<'\n'<<"next occurance "<<" after " <<NodeID[i].pattern<<(l+1)<<'\n';

//Node[count].push_back(make_Node(NodeID[i].pattern,multiplicity-l,l,0));
for(int k=0;k<j;k++)
{
//cout<<str<<'\n';
int store=find_occurrance_num(txt[k],strlen(txt[k])+1,str,(l));
//cout<<store<<'\n';
if(store<=1000)
{



char *str=txt[k]+(store+1);

//cout<<str<<'\n';
if(strlen(str)<best_length)
best_length=strlen(str);
cout<<str<<'\n';
store_str[k]=&str[0];
}


}
Node[count].push_back(make_Node(NodeID[i].pattern,NodeID[i].mul-l,best_length));
best_length=100000;
int carry_n=0,carry_l=0,l_c=0;

set<char> str=str_intersection(store_str,j);
   for (set<char>::iterator it=str.begin(); it!=str.end(); ++it)
    {




for(int k=0;k<j;k++)
{


string s_s(1,(*it));
char *str_3 = new char[1];

strcpy(str_3,s_s.c_str());
char *str_4=new char[store_str[k].size() + 1];
strcpy(str_4,store_str[k].c_str());

int *arr=find_occurrance(str_4, strlen(str_4)+1,str_3);


if((arr[1]<=mul))
{
dis=arr[0];


mul=arr[1];

}

}
string s_s(1,(*it));
cout<<*it<<'\t'<<mul<<'\t'<<dis<<'\n';
Node[count].push_back(make_Node(s_s,mul,dis));



 dis=10000000;mul=10000000;

}
for(int k=0;k<j;k++)
{

store_str[k].clear();
}

count++;

}


}

for(int i=0; i<Node.size();i++)
{

    int temp=0;
    vector<node> Node_element_temp;

for(int l=1;l<Node[i].size();l++)

   {
   for(int j=l;j<Node[i].size();j++)

   {
   if((Node[i][l].distance<Node[i][j].distance))

   {

        Node_element_temp.push_back(make_Node(Node[i][j].pattern,Node[i][j].mul,Node[i][j].distance));
        Node[i][j].pattern="";
        Node[i][j].mul=0;

        Node[i][j].distance=0;
        Node[i][j]=(make_Node(Node[i][l].pattern,Node[i][l].mul,Node[i][l].distance));
        Node[i][l]=(make_Node(Node_element_temp[0].pattern,Node_element_temp[0].mul,Node_element_temp[0].distance));
        Node_element_temp.pop_back();

   }

    }
}
   }
int count_char=0,count_item=0;
bool visited[Mul_sum];
show(Node);
for(int i=0;i<NodeID.size();i++)
{

cout<<NodeID[i].pattern<<'\t'<<NodeID[i].mul<<'\t'<<NodeID[i].distance<<'\n';

}


 for(int kl=0;kl<j;kl++)
{
string copy=&txt[kl][0];
marker1.push_back(copy.find(Node[0][2].pattern )+1);

}

 //for(int j=0; j<Node[0].size();j++){
 LCS="";
dis1=Node[0][2].distance+1;
dis=Node[0][2].distance;
 DFS(Node,Node[0][2],LCS+Node[0][2].pattern,dis,dis1,marker1);







//}
myfile.close();

return 0;
}
