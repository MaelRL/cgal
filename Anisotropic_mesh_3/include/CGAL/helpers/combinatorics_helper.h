#ifndef COMBINATORICS_HELPER_H
#define COMBINATORICS_HELPER_H

// helper for combinatorics

template<int count>
inline bool is_same_ids(int *cids, int *dids) 
{
  int t;
  for (int i = count - 1; i >= 0; i--) 
  {
    for (int j = 1; j <= i; j++) 
    {
      if (cids[j] < cids[j - 1]) 
      {
        t = cids[j];
        cids[j] = cids[j - 1];
        cids[j - 1] = t;
      }
      if (dids[j] < dids[j - 1]) 
      {
        t = dids[j];
        dids[j] = dids[j - 1];
        dids[j - 1] = t;
      }
    }
    if (cids[i] != dids[i])
      return false;
  }
  return true;
}

#endif