#ifndef RNGSTREAM_H
#define RNGSTREAM_H
#include <string>

enum ResetType {StartStream, StartSubStream, NextSubStream};

class RngStream
{
public:
  bool anti;
  RngStream (std::string name = "");
  void Reset (ResetType where);
  static void SetPackageSeed (const unsigned long seed[6]);
  void SetSeed (const unsigned long seed[6]);
  void AdvanceState (long e, long c);
  void GetState (unsigned long seed[6]) const;
  void WriteState () const;
  void WriteStateFull () const;
  double RandU01 ();
  double RandU01d ();
  long RandInt (long i, long j);

private:
  static double nextSeed[6];
  double Ig[6], Bg[6], Cg[6];
  std::string name;
};

#endif
