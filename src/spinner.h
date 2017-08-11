#pragma once

#include <iostream>
#include <iomanip>

inline void init_spinner(){
    std::cout << " ";    
}

inline void update_spinner(size_t i,size_t max){

    if (i!=0)
        std::cerr << "\b\b\b\b";            

    // Print a little spinner thing to show processing of long task:
    switch (i%4){
    case 0:{std::cerr << "\b|" ;break;}
    case 1:{std::cerr << "\b/" ;break;}
    case 2:{std::cerr << "\b-" ;break;}
    case 3:{std::cerr << "\b\\";break;}
    }

    int percentage=(i+1)*100/max;
    std::cerr << " " << std::setfill(' ') << std::setw(2) << percentage << "%";

}

inline void destroy_spinner(){
    std::cerr << "\b\b\b\b\b";
    std::cout << "\b \b";
}
