#include <fstream>
#include <iostream>

int main()
{
    std::ofstream fs("1.txt"); 

    if(!fs)
    {
        std::cerr<<"Cannot open the output file." << std::endl;
        return 1;
    }
    
    for (int i = 0; i < 10000000; i++)
    {
        fs << "1" << "\n";
    }
    fs.close();
    return 0;
}