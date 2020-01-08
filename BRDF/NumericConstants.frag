#ifndef NUMERIC_CONSTANTS_frag
#define NUMERIC_CONSTANTS_frag

#define FLT_MANT_BITS       ( 23 )
#define FLT_MIN             ( 1.175494351e-38f )
#define FLT_MAX             ( 3.402823466e+38f )
#define FLT_MAX_EXP         (  128 )
#define FLT_MIN_EXP         ( -125 )
#define FLT_EPSILON         ( 1.192092896e-07f )

#define HLF_MANT_BITS       ( 10 )
#define HLF_MIN             ( 1.0 / ( 1 << 14 ) )
#define HLF_MAX             ( 65504.0 )

#define RGBE_MAX_EXP        ( 16 )
#define RGBE_EXP_OFFSET     ( 15 )
#define RGBE_MANT_BITS      ( 9 )
#define RGBE_EXP_BITS       ( 5 )
#define RGBE_MANT_RANGE     ( 1 << RGBE_MANT_BITS )
#define RGBE_EXP_RANGE      ( 1 << RGBE_EXP_BITS )
#define RGBE_MIN_NORMALIZED ( 1.0 / ( 1 << RGBE_EXP_OFFSET ) )

#endif
