#include <iostream>

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " arg1 arg2\n";
    exit(1);
  }
  
  return 0;
}
