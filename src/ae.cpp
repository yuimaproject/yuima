#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
bool is_zero (std::string const &x){
  return((x=="0") || (x=="(0)") || (x=="{0}") || (x=="(0+0i)") || (x==""));
}



// [[Rcpp::export]]
std::vector< std::vector<std::string> > cpp_split(std::vector<std::string> const &str, std::string const sep) {
  
  int n = str.size();
  std::vector< std::vector<std::string> > out(n);
  
  for(int i=0; i<n; i++){
    
    size_t start;
    size_t end = 0;
    
    while ((start = str[i].find_first_not_of(sep, end)) != std::string::npos)
    {
      end = str[i].find(sep, start);
      out[i].push_back(str[i].substr(start, end - start));
    }
    
  }
  
  return(out);
}


// [[Rcpp::export]]
std::vector<std::string> cpp_paste(std::vector<std::string> const &x, std::vector<std::string> const &y, std::string const sep) {
  
  int n_x = x.size();
  int n_y = y.size();
  if(n_x!=n_y) Rcpp::stop("x and y must share the same length");
  
  bool is_p = (sep==" * ");
  bool is_s = (sep==" + ");
  
  std::vector<std::string> out(n_x);
  
  for(int i=0; i<n_x; i++) {
    
    out[i] = x[i] + sep + y[i];
    
    if(is_p){ 
      if(is_zero(x[i]) || is_zero(y[i])){ 
        out[i] = "0";
      }
    }
    
    if(is_s) {
      if(is_zero(y[i]) && !is_zero(x[i])) 
        out[i] = x[i];
      else if(is_zero(x[i]) && !is_zero(y[i])) 
        out[i] = y[i];
      else if(is_zero(x[i]) && is_zero(y[i]))
        out[i] = "0";
    } 
    
  }
  
  return(out);
}



// [[Rcpp::export]]
std::string cpp_collapse(std::vector<std::string> const &str, std::string const sep) {
  
  int n = str.size();
  std::string s = str[0];
  
  if(n>1){
    bool is_s = (sep==" + ");
    for(int i=1; i<n; i++) if(str[i]!="") {
      if(!is_s){
        s += sep + str[i];
      }
      else if(!is_zero(str[i])){
        if(is_zero(s)) 
          s = str[i];
        else 
          s += sep + str[i]; 
      }
    }
  }
  
  return(s);
}



// [[Rcpp::export]]
std::vector<std::string> cpp_outer(std::vector<std::string> const &x, std::vector<std::string> const &y) { 
  
  int n_x = x.size();
  int n_y = y.size();
  
  std::vector<std::string> out(n_x*n_y);
  
  int k = 0;
  for(int i=0; i<n_y; i++) for(int j=0; j<n_x; j++) {
    
    if(is_zero(x[j]) || is_zero(y[i])) 
      out[k] = "0";
    else 
      out[k] = x[j] + " * " + y[i];
    
    k++;
  }
  
  return(out);
} 



// [[Rcpp::export]]
std::vector<std::string> cpp_ito_outer(std::vector<std::string> const &x, std::vector<std::string> const &y) { 
  
  int n_x = x.size();
  int n_y = y.size();
  std::vector<std::string> out(n_x*n_y);
  
  std::vector< std::vector<std::string> > x_adds = cpp_split(x, " + ");
  std::vector< std::vector<std::string> > y_adds = cpp_split(y, " + ");
  
  int k = 0;
  for(int i=0; i<n_y; i++) for(int j=0; j<n_x; j++) {
    out[k] = cpp_collapse(cpp_outer(x_adds[j], y_adds[i]), " + ");
    k++;
  }
  
  return(out);
} 

// [[Rcpp::export]]
std::string cpp_label(std::vector<int> I){
  
  std::string s = std::to_string(I[0]);
  
  int n = I.size();
  if(n>1) for(int i=1; i<n; i++) {
    s += "," + std::to_string(I[i]);
  }
  
  return(s);
}




// [[Rcpp::export]]
std::vector<std::string> cpp_ito_product(std::vector<int> const &idx, List const &dZ, List const &Z_K, std::vector<int> const &K, int d, int a, int p, int q = 0) { 
  
  a = a-1;
  p = p-1;
  q = q-1;
  
  std::vector<std::string> tmp = std::vector<std::string> ();
  std::vector<std::string> dZ_p;
  std::vector<std::string> dZ_q;
  std::vector<int> I;
  
  // K length
  int l = K.size();
  
  // dZ.p & dZ.q
  tmp  = as< std::vector<std::string> >(dZ[ K[p]-1 ]);  
  dZ_p = std::vector<std::string>(tmp.begin() + d*a, tmp.begin() + d*(a+1));
  if(q>=0){
    tmp  = as< std::vector<std::string> >(dZ[ K[q]-1 ]);  
    dZ_q = std::vector<std::string>(tmp.begin() + d*a, tmp.begin() + d*(a+1));
  }
  
  // up to index p
  if(p==0) {
    tmp = dZ_p;
  }
  else {
    I = std::vector<int>(K.begin(), K.begin() + p);
    tmp = cpp_ito_outer(Z_K[cpp_label(I)], dZ_p);
  }
  
  // up to index q
  if(q>=0){
    if(q==(p+1)){
      tmp = cpp_ito_outer(tmp, dZ_q);
    } 
    else {
      I = std::vector<int>(K.begin() + p+1, K.begin() + q);
      tmp = cpp_ito_outer(tmp, Z_K[cpp_label(I)]);
      tmp = cpp_ito_outer(tmp, dZ_q);
    }
  }
  
  // up to index l
  int i = p;
  if(q>p) i = q;
  if(l>i+1){
    I = std::vector<int>(K.begin() + i+1, K.begin() + l);
    tmp = cpp_ito_outer(tmp, Z_K[cpp_label(I)]);
  }
  
  // return only specific indexes
  std::vector<std::string> sub;
  for(unsigned int i=0; i<idx.size(); i++) {
    sub.push_back(tmp[idx[i]]);
  }
  return(sub);
}


// [[Rcpp::export]]
std::vector<std::string> cpp_E(std::vector<std::string> str){
  
  std::vector<std::string> E;
  std::vector< std::vector<std::string> > x = cpp_split(str, " * ");
  
  int n = x.size();
  for(int i=0; i<n; i++){
    
    if (std::find(x[i].begin(), x[i].end(), "0") != x[i].end()) {
      
      E.push_back("0");
      
    } else {
      
      x[i].erase(std::remove(x[i].begin(), x[i].end(), "1"), x[i].end());
      x[i].erase(std::remove(x[i].begin(), x[i].end(), ""), x[i].end());
      
      if(x[i].size()==0) {
        E.push_back("1");
      } else {
        std::sort(x[i].begin(), x[i].end());  
        E.push_back(cpp_collapse(x[i], "_"));
      }
      
    }
    
  }
  
  return(E);
}



// [[Rcpp::export]]
List cpp_ito(List const &K_set, List const &dZ, List const &Z_K, int d, int r){
  
  std::vector<std::string> rhs;
  std::vector<std::string> lhs;
  
  for(int i=0; i<K_set.size(); i++){
    
    std::vector<std::string> tmp;
    std::vector<int> K = K_set[i];
    int l = K.size();
    
    // lhs
    tmp = cpp_E(Z_K[cpp_label(K)]);
    std::vector<int> idx;
    for(unsigned int j=0; j<tmp.size(); j++){
      if(std::find(lhs.begin(), lhs.end(), tmp[j]) == lhs.end()) {
        // store lhs
        lhs.push_back(tmp[j]);
        // store lhs idx
        idx.push_back(j);
      }
    }
    
    // rhs
    if(idx.size()>0){
      
      for(int p=1; p<=l; p++){
        
        if(p==1) tmp = cpp_ito_product(idx, dZ, Z_K, K, d, 1, p);
        else     tmp = cpp_paste(tmp, cpp_ito_product(idx, dZ, Z_K, K, d, 1, p), " + ");
        
        if(p<l) for(int q=(p+1); q<=l; q++) for(int a=2; a<=r; a++) {
          tmp = cpp_paste(tmp, cpp_ito_product(idx, dZ, Z_K, K, d, a, p, q), " + ");
        }
        
      }
      
      rhs.insert(rhs.end(), tmp.begin(), tmp.end());
      
    }
    
  }
  
  int n = rhs.size();
  std::vector< std::vector<std::string> > rhs_var(n);
  std::vector< std::vector<std::string> > adds = cpp_split(rhs, " + ");
  
  // for each equation 
  for(int i=0; i<n; i++){
    
    std::string str;
    std::vector<std::string> coef(1);
    std::vector<std::string> var(1);
    std::map<std::string,  std::map<std::string, int> > map;
    
    // for each addend 
    for(unsigned int j=0; j<adds[i].size(); j++) {
      
      int pos = 0;
      var[0]  = "";
      coef[0] = "";
      
      // parse addend
      int l = adds[i][j].length();
      if(l>0) {
        
        for(int k=0; k<l; k++) {
          
          str = adds[i][j][k];
          
          if(str=="{") {
            if(k!=0) var[0] += adds[i][j].substr(pos, k-pos);
            pos = k+1;
          }
          else if(str=="}") { 
            if(coef[0]!="") coef[0] += " * ";
            coef[0] += adds[i][j].substr(pos, k-pos);
            pos = k+1;
          }
          
        }
        
        if(str!="}") var[0] += adds[i][j].substr(pos);
        
        var = cpp_E(var);
        map[var[0]][coef[0]]++;
        
      }
      
    }
    
    
    // rhs
    rhs[i] = "";
    if(!map.empty()){
      
      std::map<std::string, std::map<std::string, int> >::iterator it1;
      std::map<std::string, int>::iterator it2;
      
      for ( it1 = map.begin(); it1 != map.end(); it1++ ) {
        if(it1->first!="0"){
          
          // simplify rhs  
          str = "";
          //std::to_string
          for ( it2 = it1->second.begin(); it2 != it1->second.end(); it2++ ) {
            if(str=="") str = std::to_string(it2->second) + " * " + it2->first;
            else str += "+" + std::to_string(it2->second) + " * " + it2->first;
          }
          if(rhs[i]=="") rhs[i] = "(" + str + ") * " + it1->first;
          else rhs[i] += " + (" + str + ") * " + it1->first;
          
          // store unique var contained in rhs
          if(it1->first!="1") rhs_var[i].push_back(it1->first);
          
        }
      }
      
    }
    
  }
  
  return(List::create(Named("lhs") = lhs , Named("rhs") = rhs, Named("rhs.var") = rhs_var));
}

