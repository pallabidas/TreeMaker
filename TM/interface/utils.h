#ifndef __UTILS_H_
#define __UTILS_H_
class PtGreater {
  public:
  template <typename T> bool operator () (const T& i, const T& j) {
    return (i.pt() > j.pt());
  }
};


class EnergyGreater {
  public:
  template <typename T> bool operator () (const T& i, const T& j) {
    return (i.energy > j.energy);
  }
};

#endif

