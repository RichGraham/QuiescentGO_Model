#include <list>

using namespace std;



struct stem{
  int top, bot;
  int Ztop, Zbot;
  int species;
  double addTop, removeTop, addBot, removeBot;
};

typedef struct stem Stem;
typedef Stem *StemPtr;
