#include <iostream>

auto errMsg =
"ERROR: The rlibmap utility is not available in ROOT6. In order to produce "
"rootmap files you can use the genreflex utils (genreflex -h for more "
"information) or the rootcling utility (rootcling -h for more information).";

int main(){
   std::cerr << errMsg << std::endl;
   return 1;
}
