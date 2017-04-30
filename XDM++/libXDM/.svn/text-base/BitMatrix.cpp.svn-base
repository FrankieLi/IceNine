#include "BitMatrix.h"

namespace GeneralLib
{

  template<  >
  BitMatrix::reference operator+=( BitMatrix::reference LHS,  const bool & RHS  )
  {
    LHS = bool( LHS ) || RHS;
    return LHS;
  }
}
