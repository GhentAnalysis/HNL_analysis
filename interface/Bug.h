#ifndef Bug_H
#define Bug_H

//include c++ library classes
#include <string>
#include <vector>
#include <iostream>

class Bug{

    public:
        Bug();
        
        void printBug( std::ostream& os = std::cout ) const;

    private:
        std::vector< std::string > lines;

};
#endif
