#include "spectrum.h"
#include "batch.h"
#include "filters.h"
#include "projectionutils.h"
#include "PeptideSpectrumMatch.h"
#include "ParameterList.h"
#include "PWizInterface.h"

#include "FilterableGraph.h"
#include "SpecnetsGraph.h"
#include "Node.h"
#include "SpectrumPairSet.h"

#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;
using namespace specnets;

template<class InputIterator> void fun1(InputIterator it, InputIterator endIt) {
  vector<float> dataV;
  //  cout<<(it==endIt)<<endl;
  
  it->serialize(dataV);
  cout<<dataV.size()<<endl;
}

void test_graph_algorithms(){
	FilterableGraph graph;// = new FilterableGraph();
	
	Node * node1 = graph.addNode();
	Node * node2 = graph.addNode();
        Node * node3 = graph.addNode();
        Node * node4 = graph.addNode();
	Edge * edge1 = graph.addEdge(node1, node2);
        Edge * edge2 = graph.addEdge(node3, node2);
	
	cout<<"Connected Components: " << graph.countComponents()<<endl;
	
	vector<Node*> node_list;
	vector<Edge*> edge_list;
	
	//int get_connected_ret = graph.getConnectedComponent(node1, node_list, edge_list);
        
        vector<int> component_sizes;
        vector<unsigned int> components_node;
        graph.getComponentSizes(component_sizes, components_node);
	
	//cout<<get_connected_ret<<endl;
	
	//cout<<"Node Size in Component: " <<node_list.size()<<endl;
	//cout<<"Edge Size in Component: " <<edge_list.size()<<endl;
        
        for(int i = 0; i < component_sizes.size(); i++){
            cout<<"component size: "<< component_sizes[i]<<endl;
        }
	
}

void test_min_component_size(char ** argv){
    cout<<"Usage: <input pairset file>"<<endl;
    
    SpectrumPairSet pairs_set;
    pairs_set.loadFromBinaryFile(argv[1]);
    
    cout<<pairs_set.size()<<endl;
    for(int i = 0; i < pairs_set.size(); i++){
        cout<<"Original Pairs: "<<pairs_set[i].spec1<<"\t"<<pairs_set[i].spec2<<endl;
    }
    
    vector<int> component_sizes;
    vector<unsigned int> components_node;
    SpecnetsGraph spec_graph(pairs_set);
    
    spec_graph.getComponentSizes(component_sizes, components_node);
    
    for(int i = 0; i < component_sizes.size(); i++){
        cout<<"component size: "<< component_sizes[i]<<endl;
    }
    
    spec_graph.filter_graph_component_size(3);
    cout<<"FILTERING COMPONENTS"<<endl;
    component_sizes.clear();
    components_node.clear();
    spec_graph.getComponentSizes(component_sizes, components_node);
    
    for(int i = 0; i < component_sizes.size(); i++){
        cout<<"component size: "<< component_sizes[i]<<endl;
    }
    
    vector<unsigned int> deleted_edges = spec_graph.get_pairs_deleted();
    
    for(int i = 0; i < deleted_edges.size(); i++){
        cout<<"DELETED EDGE: "<<deleted_edges[i]<<endl;
    }
}

void write_test_output_spectrumpairsset(char ** argv){
    cout<<"Usage: <output pairset file>"<<endl;
    
    SpectrumPairSet pairs_set;
    SpectrumPair pair1;
    pair1.spec1 = 1;
    pair1.spec2 = 2;
    pair1.score1 = 0.7;
    
    SpectrumPair pair2;
    pair2.spec1 = 2;
    pair2.spec2 = 3;
    pair2.score1 = 0.7;
    
    SpectrumPair pair3;
    pair3.spec1 = 3;
    pair3.spec2 = 4;
    pair3.score1 = 0.5;
    
    pairs_set.push_back(pair1);
    pairs_set.push_back(pair2);
    pairs_set.push_back(pair3);
    
    pairs_set.saveToBinaryFile(argv[1]);
}

void pwiz_to_mgf(char ** argv){
    cout<<"Usage: <input file> <output mgf>"<<endl;
    SpecSet spec1;
    
    PWizInterface pwiz;
    int loadOk = pwiz.loadDataUsingPWiz(argv[1], spec1, 2);
    spec1.SaveSpecSet_mgf(argv[2]);
}

void test_mgf_loadsave(char ** argv){
    cout<<"Usage: <input file> <output mgf>"<<endl;
    SpecSet spec1;
    spec1.Load(argv[1]);
    
    spec1.removeSpectraBelowMinPeaks(5);
    
    spec1.SaveSpecSet_mgf(argv[2]);
}

void test_mgf_loadscans(char ** argv){
    cout<<"Usage: <input file> <output mgf> <scan>"<<endl;
    SpecSet spec1;
    cout<<"LOADING"<<endl;
    
    spec1.LoadSpecSet_mgf(argv[1], atoi(argv[3]), -1);
    spec1.SaveSpecSet_mgf(argv[2]);
}

void fix_mgf_scans(char ** argv){
    cout<<"Usage: <input file> <output mgf>"<<endl;
    SpecSet spec1;
    spec1.Load(argv[1]);
    
    for(int i = 0; i < spec1.size(); i++){
        spec1[i].scan = i + 1;
    }
    
    spec1.SaveSpecSet_mgf(argv[2]);
}

void test_outputPM(char ** argv){
    cout<<"Usage: <input file>"<<endl;
    SpecSet spec1;
    spec1.Load(argv[1]);
    cout<<spec1.size()<<endl;
    for(int i = 0; i < spec1.size(); i++){
        cout<<i<<"\t"<<spec1[i].scan<<"\t"<<spec1[i].parentMass<<endl;
    }
}

void test_networking_cosine(char ** argv){
    cout<<"Usage: <input mgf> <scan1> <scan2> <peakTol>"<<endl;
    SpecSet spec1;
    spec1.Load(argv[1]);
    
    
    
    unsigned int matched_peaks;
    float peakTol = atof(argv[4]);
    float cosine = spec1[atoi(argv[2])-1].scoreMatch(spec1[atoi(argv[3])-1], peakTol, matched_peaks, false);
    //spec1[atoi(argv[2])-1].output(cout);
    
    cout<<"Cosine:\t" <<cosine<<endl<<"Matched Peaks:\t"<<matched_peaks<<endl;
    
}

void create_organize_scans(char ** argv){
    cout<<"Usage: <input mgf> <output mgf>"<<endl;
    SpecSet spec1;
    spec1.Load(argv[1]);
    cout<<spec1.size()<<endl;
    spec1.SaveSpecSet_mgf(argv[2]);
}

void create_FDA_library(char ** argv){
    cout<<"Usage: <input mgf> <output mgf> <library name>"<<endl;
    SpecSet spec1;
    spec1.Load(argv[1]);
    
    cout<<"Loaded file of size "<<spec1.size()<<endl;
    
    string file_name(argv[1]);

    for(int i = 0; i < spec1.size(); i++){
        float retention_time = spec1[i].retention_time;
        spec1[i].instrument_name = "qtof";
        
        psmPtr psm(new PeptideSpectrumMatch);
        psm->m_annotation = "*..*";
        psm->m_ionmode = "positive";
        //psm->m_organism.push_back(argv[3]);
        //psm->m_submission_metadata.push_back(argv[1]);
        //psm->m_compound_name.push_back(strip_extension(file_name));
        //psm->m_smiles.push_back(argv[3]);
        //psm->m_InChI.push_back(argv[3]);
        //psm->m_InChI_Aux.push_back(argv[3]);
        //psm->m_notes = argv[1];
        
        spec1[i].psmList.push_back(psm);
    }


    spec1.SaveSpecSet_mgf(argv[2]);
}

void create_unannotated_library(char **argv){
    //Ming's test library conversion code
    cout<<"Usage: <input pklbin> <output mgf> <organism name>"<<endl;
    SpecSet spec1;
    spec1.Load(argv[1]);
    
    cout<<"Loaded file of size "<<spec1.size()<<endl;

    for(int i = 0; i < spec1.size(); i++){
        float retention_time = spec1[i].retention_time;
        float mz = spec1[i].parentMZ;
        int scan = spec1[i].scan;
        
        
        
        
        psmPtr psm(new PeptideSpectrumMatch);
//         psm->m_organism.push_back(argv[3]);
//         psm->m_submission_metadata.push_back(argv[1]);
//         psm->m_compound_name.push_back(argv[3]);
//         psm->m_smiles.push_back(argv[3]);
//         psm->m_InChI.push_back(argv[3]);
//         psm->m_InChI_Aux.push_back(argv[3]);
//         psm->m_notes = argv[1];
//         
        spec1[i].psmList.push_back(psm);
        spec1[i].ITOL = retention_time;
    }


    spec1.SaveSpecSet_mgf(argv[2]);
}

int main(int argc, char **argv) {
  //SpecSet specs;

  //specs.loadPklBin(argv[1]);
  //specs.SaveSpecSet_mgf(argv[2]);

//   vector<Results_PA> pairsPA;
//   vector<Results_ASP> pairsASP;
//   int idx;
// 
//   cout<<"Load_resultsASPbin: "<<Load_resultsASPbin("testASP.bin",pairsASP)<<endl;
//   cout<<"  -> "<<pairsASP.size()<<" pairs";
//   idx = (int)pairsASP.size()-1;
//   if(idx<0) cout<<endl; else cout<<", last pair = ("<<pairsASP[idx].spec1<<","<<pairsASP[idx].spec2<<","<<pairsASP[idx].shift1<<","<<pairsASP[idx].score1<<","<<pairsASP[idx].score2<<")\n";
//   Save_resultsASPbin("testASPout1.bin",pairsASP);
//   
//   cout<<"Load_results_bin: "<<Load_results_bin("testASP.bin",pairsASP)<<endl;
//   cout<<"  -> "<<pairsASP.size()<<" pairs";
//   idx = (int)pairsASP.size()-1;
//   if(idx<0) cout<<endl; else cout<<", last pair = ("<<pairsASP[idx].spec1<<","<<pairsASP[idx].spec2<<","<<pairsASP[idx].shift1<<","<<pairsASP[idx].score1<<","<<pairsASP[idx].score2<<")\n";
//   Save_results_bin("testASPout2.bin",pairsASP.size(),pairsASP.begin());*/
//   
  //Ming's Test Cosine Code
  
  //fix_mgf_scans(argv);
  //pwiz_to_mgf(argv);
  //test_outputPM(argv);
  //pwiz_to_mgf(argv);
  //test_mgf_loadsave(argv);
  //test_mgf_loadscans(argv);
  //test_mgf_loadsave(argv);
  //create_organize_scans(argv);
  //test_graph_algorithms();
  write_test_output_spectrumpairsset(argv);
  test_min_component_size(argv);
  return 0;
  
//   SpecSet spec1;
//   SpecSet spec2;
//   
//   spec1.Load(argv[1]);
//   spec2.Load(argv[3]);
//   
//   float cosine = full_spectrum_similarity(spec1[atoi(argv[2])+1], spec2[atoi(argv[4])+1]);
// 
//   cout<<"COSINE: " <<cosine<<endl;*/
  
  
  //create_FDA_library(argv);
  //create_organize_scans(argv);
  //create_unannotated_library(argv);
  //test_networking_cosine(argv);
  //test_outputPM(argv);
  
  /*
  //Ming's Test for bipartite cosine in spectrum
  cout<<"Usage <input mgf> <index 1> <index 2>"<<endl;
  SpecSet spec1;
  spec1.Load(argv[1]);
  int index1 = atoi(argv[2]);
  int index2 = atoi(argv[3]);
  unsigned int matchedpeaks;
  
  float score = spec1[index1].scoreMatch(spec1[index2], 0.05, matchedpeaks);
  cout<<score<<"\t"<<matchedpeaks<<endl;*/
  
  
  //   cout<<"COSINE: " <<cosine<<endl;
  
}
