#include <vector>           // Pouzivej knihovnu pro pole
 
using namespace std;        // Hledej objekty v knihovne std (std::vector -> vector)


int main(int argc, char** argv) {

  vector<double> p(100);    // pole o delce 100 prvku
  
  for (int i=0; i<100; i++) {
    p[i] = 0.5*i + 1;
  }

  for (int i=0; i<p.size(); i++) {   // p.size() vraci delku pole
    p[i] = 0.5*i + 1;
  }

  return 0;
}
