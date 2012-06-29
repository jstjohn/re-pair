/*
  Modified source from:
  http://code.google.com/p/ngopt/source/browse/trunk/tools/pair_reads/repair.cpp?r=85
  "repair"
*/
#include <istream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <cstdio>
#include <cstring>

using namespace std;

struct read {
  string hdr;
  string comment;
  string seq;
  string qual;
  int pair;
};


int pair_reads_call = 0;
map<string,struct read> pair1;
map<string,struct read> pair2;

void pair_reads(istream& in, bool fastq, bool casava){ 
  pair_reads_call++;
  struct read r;
  string hdr;
  string q_hdr;
  string seq;
  string qual;
  string key;
  string comment;
  string line1;
  string line2;
  string line3;
  string line4;
  while(getline(in, line1) && getline(in, line2) && getline(in, line3) && getline(in, line4)){
    stringstream lss1(line1);
    lss1 >> hdr;
    if(casava)
       lss1 >> comment;
    stringstream lss2(line2);
    lss2 >> seq;
    if (fastq) {
      stringstream lss3(line3);
      lss3 >> q_hdr;
      stringstream lss4(line4);
      lss4 >> qual;
    }
    hdr = hdr.substr(1); //leave off the initial @
    if(casava)
      r.hdr = hdr + "\t" + comment;
    else
      r.hdr = hdr;
    r.seq = seq;

    if (fastq){
      r.qual = qual;
    }

    bool read1 = false;
    bool read2 = false;

    if(casava){
      //casava style fastq headers:
      //@HWI-ST593:192:D0D55ABXX:8:1101:1145:2077 1:N:0:
      
      key = hdr;
      if(comment.at(0) == '1')
	read1 = true;
      if(comment.at(0) == '2' || comment.at(0) == '3')
	read2 = true;

    }else{ // standard fastq headers (/1 and /2 (or 3)) with no comment
      key = hdr.substr(0,hdr.length()-2);
      if (hdr.at(hdr.length()-1) == '1')
	read1 = true;
      if (hdr.at(hdr.length()-1) == '2' || hdr.at(hdr.length()-1) == '3')
	read2 = true;

    }

    if (read1) {
      r.pair = 1;
      pair1[key] = r;
    } else if (read2) {
      r.pair = 2;
      pair2[key] = r;
    } else {
      cerr << "Unable to pair read: >>" << r.hdr << "<<" << endl;
    }
  }
}

void print_paired(string prefix, string base, string suffix, bool fastq){

  ofstream p1out((prefix+base+"_p1"+suffix).c_str());
  ofstream p2out((prefix+base+"_p2"+suffix).c_str());
  ofstream upout;

  map<string, struct read>::iterator it;
  struct read*  tmp1;
  struct read*  tmp2;
  for (it=pair1.begin(); it!=pair1.end(); it++) {
    if (pair2.find(it->first) != pair2.end()) {
      tmp1 = & it->second;
      tmp2 = & pair2.find(it->first)->second; 
      if (fastq) {
	p1out << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
	p2out << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
      } else {
	p1out << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
	p2out << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;

      }
      pair2.erase(it->first);
    } else {
      if (!upout.is_open()){
	upout.open((prefix+base+"_up"+suffix).c_str());
      }
      tmp1 = & it->second;
      if (fastq)
	upout << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
      else
	upout << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
    }
  }

  for (it=pair2.begin(); it!=pair2.end(); it++) {
    if (!upout.is_open()){
      upout.open((prefix+base+"_up"+suffix).c_str());
    }
    tmp2 = & it->second;
    if (fastq)
      upout << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
    else 
      upout << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;
  }

  p1out.close();
  p2out.close();
  if (upout.is_open())
    upout.close();

}

void print_shuffled(string prefix, string base, string suffix, bool fastq){

  ofstream paired((prefix+base+"_shuf"+suffix).c_str());
  ofstream unpaired;

  map<string, struct read>::iterator it;
  struct read*  tmp1;
  struct read*  tmp2;
  for (it=pair1.begin(); it!=pair1.end(); it++) {
    if (pair2.find(it->first) != pair2.end()) {
      tmp1 = & it->second;
      tmp2 = & pair2.find(it->first)->second; 
      if (fastq) {
	paired << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
	paired << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
      } else {
	paired << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
	paired << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;

      }
      pair2.erase(it->first);
    } else {
      if (!unpaired.is_open()) {
	unpaired.open((prefix+base+"_up"+suffix).c_str());
      }
      tmp1 = & it->second;
      if (fastq)
	unpaired << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
      else
	unpaired << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
    }
  }

  for (it=pair2.begin(); it!=pair2.end(); it++) {
    if (!unpaired.is_open()) {
      unpaired.open((prefix+base+"_up"+suffix).c_str());
    }
    tmp2 = & it->second;
    if (fastq)
      unpaired << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
    else 
      unpaired << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;
  }

  paired.close();
  if (unpaired.is_open())
    unpaired.close();

}

void pipe_seq(istream& in, ostream& out, bool fastq) {
  string hdr;
  string seq;
  if (in >> hdr){
    in >> seq;
    out << hdr << '\n' << seq << endl;
    if (fastq) {
      in >> hdr;
      in >> seq;
      out << hdr << '\n' << seq << endl;
    }
  }
}

void shuffle_paired(ifstream& in1, ifstream& in2, ofstream& out, bool fastq) {
  while (in1.good() && in2.good()){
    pipe_seq(in1,out,fastq);
    pipe_seq(in2,out,fastq);
  }
}

void split_shuffled(istream& in, ofstream& p1out, ofstream& p2out, bool fastq){
  while (in.good()){
    pipe_seq(in,p1out,fastq);
    pipe_seq(in,p2out,fastq);
  }
}

void usage(char* name){
  cout << "Usage: "<<name<<" [options] <base> <read1.fq> <read2.fq> ... <readN.fq[.gz]>"<<endl;
  cout << " where "<< endl;
  cout << "       <base>       basename for output files\n";
  cout << "       <reads1..N.fq>   The fastq files to repair, if none provided\n";
  cout << "                     then standard in is used.\n";
  cout << " options:\n";
  cout << "        --casava    input files are in illumina casava v1.8 fastq format.\n";
  cout << "        -p <string> a prefix to add to <base>. This can be used to\n"; 
  cout << "                    specify an output directory.\n";
  cout << "        -s <string> the suffix to append to the output files.\n"; 
  cout << "        --shuf      print pairs in one file where paired reads are printed\n";
  cout << "                    on consecutive lines (a.k.a. shuffled).\n";
  cout << "        --paired    assume reads are paired.\n";
  cout << "        --quiet     do not print progress messages.\n";
  cout << "        --debug     run in debug mode.\n\n";

}

int main(int argc, char** argv){
  if (argc == 1) {
    usage(argv[0]);
    return 0;
  }
  string prefix = "";
  string suffix = "";
  string base = "";
  bool fastq = true;
  bool shuffle = false;
  bool paired = false;
  bool debug = false;
  bool casava = false;
  bool quiet = false;
  bool gzip = false;
  int start = 1;
  int i = 1;
  while (argv[i][0] == '-') {
    if (argv[i][1]=='s'){
      suffix = argv[++i];
      start+=2;
    } else if (argv[i][1]=='p') {
      prefix = argv[++i];
      start+=2;
    } else if (argv[i][1]=='-'){
      if (strcmp(argv[i],"--shuf")==0){
	shuffle = true;
      } else if (strcmp(argv[i],"--paired")==0){
	paired = true;
      } else if (strcmp(argv[i],"--casava")==0){ 
	casava = true;
      }else if (strcmp(argv[i],"--quiet")==0){
	quiet = true;
      } else if (strcmp(argv[i],"--debug")==0){
	debug = true;
      } else if (strcmp(argv[i],"--gzip")==0){
	gzip = true;
      }
      start++;
    } else {
      cerr << "Unrecognized argument: " << argv[i] << endl;
    }
    i++;
  }
  if (debug) cerr << argc << " arguments. Starting at " << start << endl;

  char c;
  istream* in;
  base = argv[start++];   
  if (fopen(base.c_str(),"r")) {
    cerr << "Missing <base> argument\n";
    usage(argv[0]);
    return 0;
  }
  if (paired) { 
    if (shuffle) { // need two files
      if (argc - start == 2) {
	ifstream in1(argv[start++]);
	ifstream in2(argv[start++]);
	ofstream out((prefix+base+"_shuf"+suffix).c_str());
	fastq = in1.peek()=='@';
	if (!quiet) cout << "Shuffling paired reads from " << argv[start-2] << " and " << argv[start-1] << "\n";
	shuffle_paired(in1,in2,out,fastq);
      } else {
	cerr << "Two files are required for shuffling paired reads.\n";
	usage(argv[0]);
      }
    } else {  
      if (!quiet) cout << "Splitting shuffled reads from";
      ofstream p1out((prefix+base+"_p1"+suffix).c_str());
      ofstream p2out((prefix+base+"_p2"+suffix).c_str());
      if (argc - start == 0) {
	in = &cin;
	fastq = in->peek() == '@';
	if (!quiet) cout << " standard input.\n";
	split_shuffled(*in,p1out,p2out,fastq);
      } else {
	filebuf fb;
	do {                    
	  fb.open(argv[start++],ios::in);
	  if (!quiet) cout << "\n\t" << argv[start-1];
	  in = new istream(&fb);
	  fastq = in->peek() == '@';
	  //      if (debug) cerr << start << "  " << argv[start] << endl;
	  split_shuffled(*in,p1out,p2out,fastq);
	  delete in;
	} while (start < argc);
	if (!quiet) cout << '\n';
	//      fb.open(argv[i],ios::in);
	//      in = new istream(&fb);
	//      split_shuffled(*in,p1out,p2out);

      }
      return 1;
    } 
  } else  { // reads need to be re-paired using a hash
    if (!quiet) cout << "Pairing reads from" ;
    if (argc - start == 0) {
      in = &cin;
      fastq = in->peek() == '@';
      if (!quiet) cout << " standard input.\n";
      pair_reads(*in,fastq,casava);
    } else {
      filebuf fb;
      do {
	fb.open(argv[start++],ios::in);
	if (!quiet) cout << "\n\t" << argv[start-1];
	in = new istream(&fb);
	fastq = in->peek() == '@';
	pair_reads(*in,fastq,casava);
	fb.close();
	delete in;
      } while (start < argc);
      if (!quiet) cout << '\n';
      //                      for (int i = start; i < argc; i++){
      //                              cerr << "Loading reads from " << argv[i] << endl;
      //                              fb.open(argv[i],ios::in);
      //                              in = new istream(&fb);
      //                              pair_reads(*in);
      //                              fb.close();
      //                              delete in;
      //                      }
    }
    if (shuffle)
      print_shuffled(prefix,base,suffix,fastq);
    else
      print_paired(prefix,base,suffix,fastq);
  }

}
