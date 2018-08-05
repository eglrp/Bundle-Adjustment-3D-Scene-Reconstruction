#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "stdio.h"
#include "stdlib.h"
#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include "assert.h"
using namespace std;

struct Vertex{
    vector<int> nbrs;    
};

// v1<v2 must be guaranteed
struct Edge{
    int v1;
    int v2;
};    

int n_cams;
vector<Vertex> vertices;
vector<Vertex> Skeletal_vertices;
vector<int> discovered_vertices;
vector<Edge> edges;
vector<Edge> saved_edges;
int unassigned_vtx_count;

bool checkIsForest(const vector<Vertex> &G, int num_components)
{
    int edge_count=0;
    for(int i=0;i<G.size();i++)
        edge_count+=G[i].nbrs.size();
        
    edge_count/=2;
    cout<<"edge_count="<<edge_count<<"\n";
    return (G.size()==edge_count+num_components);
}
/*
void process_vertex(int vid)
{
    Vertex& v=vertices[vid];
    vector<int> cur_vertices;
    for(int i=0;i<v.nbrs.size();i++)
    {
        int nbr_id=v.nbrs[i];
        if(discovered_vertices[nbr_id]==0)
        {
            Edge e;
            if(vid<nbr_id){
                e.v1=vid;e.v2=nbr_id;
            }else{
                e.v1=nbr_id;e.v2=vid;
            }
            saved_edges.push_back(e);
            Skeletal_vertices[vid].nbrs.push_back(nbr_id);
            Skeletal_vertices[nbr_id].nbrs.push_back(vid);
            cur_vertices.push_back(nbr_id);
            discovered_vertices[nbr_id]=1;
            unassigned_vtx_count--;
        }    
    }
    int children_max_degree=-1;
    int cand_vtx=-1;
    for(int i=0;i<cur_vertices.size();i++)
    {
        Vertex& next_vtx=vertices[cur_vertices[i]];
        int degree=0;
        for(int j=0;j<next_vtx.nbrs.size();j++)
        {
            if(discovered_vertices[next_vtx.nbrs[j]]==0)
                degree++;
        }
        if(children_max_degree<degree)
        {
            children_max_degree=degree;
            cand_vtx=cur_vertices[i];
        }
    }
    if(cand_vtx>0)
    {
        process_vertex(cand_vtx);
        discovered_vertices[vid]=2;
    }    
}
*/
void expand_vtx(int vid, set<int> &cand_vertices, vector<int> &unmarked_deg)
{
    discovered_vertices[vid]=1;
    unassigned_vtx_count--;
    for(int j=0;j<vertices[vid].nbrs.size();j++)
    {
        int nbr_idx=vertices[vid].nbrs[j];
        if(discovered_vertices[nbr_idx]==0)
            unmarked_deg[nbr_idx]--;
    }    
    vector<int> &nbrs=vertices[vid].nbrs;
    for(int i=0;i<nbrs.size();i++)
    {
        int nbr_id=nbrs[i];
        if(discovered_vertices[nbr_id]==0)
        {
            Edge e;
            if(vid<nbr_id){
                e.v1=vid;e.v2=nbr_id;
            }else{
                e.v1=nbr_id;e.v2=vid;
            }    
            saved_edges.push_back(e);
            Skeletal_vertices[vid].nbrs.push_back(nbr_id);
            Skeletal_vertices[nbr_id].nbrs.push_back(vid);
            discovered_vertices[nbr_id]=1;
            unassigned_vtx_count--;
            cand_vertices.insert(nbr_id);
            for(int j=0;j<vertices[nbr_id].nbrs.size();j++)
            {
                int nbr_idx=vertices[nbr_id].nbrs[j];
                if(discovered_vertices[nbr_idx]==0)
                    unmarked_deg[nbr_idx]--;
            }    
        }
    }
    discovered_vertices[vid]=2;
}

void build_single_MLST()
{
    int root_vtx=-1;
    int max_degree=-1;
    for(int i=0;i<vertices.size();i++)
    {
        int deg=vertices[i].nbrs.size();    
        if((discovered_vertices[i]==0)&&(deg>max_degree))
        {
            max_degree=deg;
            root_vtx=i;
        }
    }  
    
    
    if(root_vtx>=0)
    {
        set<int> cand_vertices;
        vector<int> unmarked_deg;
        unmarked_deg.resize(vertices.size());
        for(int i=0;i<unmarked_deg.size();i++)            
            unmarked_deg[i]=vertices[i].nbrs.size();
        
        cand_vertices.insert(root_vtx);
        while(cand_vertices.size()>0)
        {
            int max_unmarked_deg=-1;
            set<int>::iterator next_vtx;
            for(set<int>::iterator it=cand_vertices.begin();it!=cand_vertices.end();it++)
            {
                int vtx_idx=*it;                       
                if(unmarked_deg[vtx_idx]>max_unmarked_deg)
                {
                    max_unmarked_deg=unmarked_deg[vtx_idx];
                    next_vtx=it;
                }                  
                expand_vtx(*it, cand_vertices, unmarked_deg);
                cand_vertices.erase(it);
            } 
        }
    }
}

void build_MLST()
{   
    Skeletal_vertices.resize(vertices.size());
    discovered_vertices.resize(vertices.size());
    unassigned_vtx_count=vertices.size();
    for(int i=0;i<vertices.size();i++)
        discovered_vertices[i]=0;
        
    saved_edges.clear();
    cout<<"Start building MLST forest\n";
    int component_count=0;
    while(unassigned_vtx_count>0)
    {
        build_single_MLST();
        component_count++;
    }       
    cout<<"Finish building MLST forest. There're a total of "<<component_count<<" components\n";
    cout<<"Is valid MLST forest?"<<checkIsForest(Skeletal_vertices, component_count)<<"\n";
}

int BFS_dist(int source, int dest)
{    
    vector<int> BFS_dists(vertices.size());
    for(int i=0;i<vertices.size();i++)
        BFS_dists[i]=-1;    
    BFS_dists[source]=0;
    
    queue<int> BFS_queue;
    BFS_queue.push(source);
    while(BFS_queue.size()>0)
    {
        int next_vtx=BFS_queue.front();
        BFS_queue.pop();
        vector<int> &nbrs=Skeletal_vertices[next_vtx].nbrs;
        for(int i=0;i<nbrs.size();i++)
        {
            int nbr_vtx=nbrs[i];            
            if(BFS_dists[nbr_vtx]==-1)
            {
                BFS_queue.push(nbr_vtx);
                BFS_dists[nbr_vtx]=BFS_dists[next_vtx]+1;
            }
            if(nbr_vtx==dest)
                return BFS_dists[dest];
        }
    }     
    return 999999;
}

inline bool edgesort_comparison_func(Edge e1, Edge e2)
{
    int score1=discovered_vertices[e1.v1]+discovered_vertices[e1.v2];
    int score2=discovered_vertices[e2.v1]+discovered_vertices[e2.v2];
    return (score1>score2);
}

void build_skeletal(int spanner_t)
{
    sort(edges.begin(), edges.end(), edgesort_comparison_func);
    for(int i=0;i<edges.size();i++)
    {
        Edge e=edges[i];
        int dist=BFS_dist(e.v1, e.v2);
        if(dist>=spanner_t)
        {             
            saved_edges.push_back(e);
            Skeletal_vertices[e.v1].nbrs.push_back(e.v2);
            Skeletal_vertices[e.v2].nbrs.push_back(e.v1);
        }
    }
}

void save_result(const char* output_filename)
{
    ofstream output;
    output.open(output_filename);
    for(int i=0;i<saved_edges.size();i++)
    {
        Edge& e=saved_edges[i];
        output<<e.v1<<" "<<e.v2-n_cams<<endl;
    }
}

void BFS(int source, vector<bool> &discovered)
{       
    queue<int> BFS_queue;
    BFS_queue.push(source);
    discovered[source]=true;
    while(BFS_queue.size()>0)
    {
        int next_vtx=BFS_queue.front();
        BFS_queue.pop();
        vector<int> &nbrs=vertices[next_vtx].nbrs;
        for(int i=0;i<nbrs.size();i++)
        {
            int nbr_vtx=nbrs[i];   
            if(!discovered[nbr_vtx])
            {
                discovered[nbr_vtx]=true;
                BFS_queue.push(nbr_vtx);
            }                     
        }
    }         
}

void init(const char* filename)
{
    ifstream input;
    input.open(filename);
    if(input.fail())
    {
        cout<<"Error opening the input file\n";
        exit(0);
    }else{
        cout<<"Successfully open input file\n";
    }    
    int n_pts, n_obs;
    input>>n_cams;
    input>>n_pts;
    input>>n_obs;
    int cid, pid;
    double obs1, obs2;
    cout<<"n_cams="<<n_cams<<" n_pts="<<n_pts<<" n_obs="<<n_obs<<"\n";
    
    vertices.resize(n_cams+n_pts);
    edges.resize(n_obs);    
    for(int i=0;i<n_obs;i++)
    {
        input>>cid;
        input>>pid;
        input>>obs1;
        input>>obs2;        
        vertices[cid].nbrs.push_back(pid+n_cams);
        vertices[pid+n_cams].nbrs.push_back(cid);     
        Edge e;e.v1=cid;e.v2=pid+n_cams;
        edges[i]=e;   
    }
    
    for(int i=0;i<vertices.size();i++)
    {
        if(vertices[i].nbrs.size()==0)
        {
            cout<<"vertex "<<i<<" has no neighbors\n";
        }
    }      
    
    vector<bool> discovered;
    discovered.resize(vertices.size());
    for(int i=0;i<vertices.size();i++)
        discovered[i]=false;
    
    int component_count=0;    
    for(int i=0;i<vertices.size();i++)
    {
        if(!discovered[i])
        {
            BFS(i, discovered);
            component_count++;
        }    
    }    
    cout<<"Input graph has "<<component_count<<" components\n";           
    input.close();
}

int main(int argc, char** argv)
{
    init(argv[1]);
    cout<<"Finish init\n";
    build_MLST();
    cout<<"Finish building MLST forest\n";
    int spanner_t=atoi(argv[3]);
    cout<<"Spanner_t="<<spanner_t<<"\n";
    build_skeletal(spanner_t);
    cout<<"Finish building skeletal graph\n";
    save_result(argv[2]);
    return 0;
}
