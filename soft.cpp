#include <iostream>

#include <climits>
#include <cfenv> //fesetround
#include <cmath> //nextafter
#include <cfloat> 

#define PRINT_FLOAT_OWN(x) std::cout << *reinterpret_cast<float*>(& x) << '\n';
#define PRINT_FLOAT_STD(x) std::cout << x << '\n';
#define PRINT_FLOAT_BIT(x) { float r = x; std::cout << *reinterpret_cast<uint32_t*>(& r) << '\n'; }
#define COMPARE_FLOAT(x, y) (*reinterpret_cast<float*>(&x) ==( y ))
#define CHECK(x, y) \
	res = x;\
	if(COMPARE_FLOAT(res, ( y )))\
		std::cout << " V";\
	else\
	{\
		std::cout << " X\n";\
		PRINT_FLOAT_OWN(res);\
		PRINT_FLOAT_STD(y);\
		std::cout << res << '\n';\
		float t = y;\
		PRINT_FLOAT_BIT(t);\
		std::cout << '\n';\
	}

#define CHECK_INT(x, y) \
	if((x) == (y))\
		std::cout << " V";\
	else\
	{\
		std::cout << " X\n";\
		std::cout << (x) << '\n';\
		std::cout << (y) << '\n';\
		std::cout << '\n';\
	}

#define IF_A_B_BOTH(a, b, x, y, z)\
	if(a)\
	{\
		if(b)\
			return z;\
		return x;\
	}\
	if(b)\
	{\
		return y;\
	}

#define QNAN 0xFFC00000
#define INF  0x7F800000
#define NINF 0xFF800000

#define BIAS  1024

/*
All references to the standard are to IEEE754 2008 

all functions are standard compliant with the exception of signaling errors

not implemented parts of standard:
	different roundings                              ; the only one available is TOWARD_ZERO
	remainder(x, y)                                  ;
	convert to/from decimal character sequence       ;
	anything related to error handling and signaling ;

	fused multiply add                               ; I am too tired to do it in the nearest future

no reason to implement 
	radix                                            ; everything is in base 2 
	is canonical                                     ; every number accessible to user thru functions other than pack/unpack is canonical
	minmag, maxmag, totalordermag                    ; not needed since user can use usual functions with abs() on arguments
	copysign(x)                                      ; x = x 
*/


/*
	internal representation is binary float with extended precision
	this way there is no need to handle subnormals in different way than normals, pack() and unpack() take care of that 
*/

enum class cls{
	snan,
	qnan,
	ninf,
	nnrm,
	nsbn,
	nzro,
	pzro,
	psbn,
	pnrm,
	pinf
};

struct Content{
	uint32_t mantisa;
	uint16_t exponent;
	uint8_t  sign;
};

/*
	operations used for internal calculation
*/
constexpr Content  unpack              (const uint32_t a);
constexpr uint32_t pack                (const Content& c);

/*
	general operation
	section 5.3.1
*/
constexpr uint32_t whl                 (const uint32_t a);
constexpr uint32_t next_Up             (const uint32_t a);
constexpr uint32_t next_Down           (const uint32_t a);
constexpr uint32_t min_Val             (const uint32_t a, const uint32_t b);
constexpr uint32_t max_Val             (const uint32_t a, const uint32_t b);

/*
	logBFormat operations
	section 5.3.3
*/
constexpr uint32_t scl                 (const uint32_t a, const int32_t  N);
constexpr int32_t  log                 (const uint32_t a);

/*
	arithmetic operations
	section 5.4.1
*/
constexpr uint32_t add                 (const uint32_t a, const uint32_t b);
constexpr uint32_t sub                 (const uint32_t a, const uint32_t b);
constexpr uint32_t mul                 (const uint32_t a, const uint32_t b);
constexpr uint32_t div                 (const uint32_t a, const uint32_t b);
constexpr uint32_t sqt                 (const uint32_t a);
constexpr uint32_t int_To_Fl           (const int32_t  i);
constexpr int32_t  fl_To_Int           (const uint32_t a);


/*
	sing bit operations
	section 5.5.1
*/
constexpr uint32_t neg                 (const uint32_t a);
constexpr uint32_t abs                 (const uint32_t a);
constexpr uint32_t cps                 (const uint32_t a, const uint32_t b);


/*
	comparisons
	section 5.6.1
	prefix q means quiet
*/
constexpr bool     qeq                 (const uint32_t a, const uint32_t b);
constexpr bool     qne                 (const uint32_t a, const uint32_t b);
constexpr bool     qgr                 (const uint32_t a, const uint32_t b);
constexpr bool     qge                 (const uint32_t a, const uint32_t b);
constexpr bool     qls                 (const uint32_t a, const uint32_t b);
constexpr bool     qle                 (const uint32_t a, const uint32_t b);
constexpr bool     total_Order         (const uint32_t a, const uint32_t b);

/*
	conformance predicates
	section 5.7.1
*/
constexpr bool     is_754_Version_1985 ();
constexpr bool     is_754_Version_2008 ();

/*
	general operations
	section 5.7.2
*/
constexpr cls      float_Class         (const uint32_t a);
constexpr bool     is_Neg              (const uint32_t a);
constexpr bool     is_Nrm              (const uint32_t a);
constexpr bool     is_Fin              (const uint32_t a);
constexpr bool     is_Zero             (const uint32_t a);
constexpr bool     is_Sbn              (const uint32_t a);
constexpr bool     is_Inf              (const uint32_t a);
constexpr bool     is_NaN              (const uint32_t a);
constexpr bool     is_sNaN             (const uint32_t a);




constexpr bool is_754_Version_1985() { return false; } //no idea, playing safe
constexpr bool is_754_Version_2008() { return false; } //reasons at the beggining of file

/*
3,3 sets of floating point data
table 3.2 parameters defining basic format floating-point numbers

this code implements binary32 format
*/
constexpr Content unpack(const uint32_t a)
{
	uint8_t sign      = a >> 31;                         //1
	uint16_t exponent = ((a >> 23) & 0xFF) + 1024 - 127;// + (t == 0); //8
	uint32_t mantisa  = (a & 0x7FFFFF); //23

	//if exponent is zero, then return actual zero
	//else return normal number with extended precision
	if(exponent == 1024 - 127)
	{
		if(mantisa == 0)
		{
			exponent = 0;
			goto ret;
		}

		mantisa <<= 1;
		//change subnormals to normals
		while(mantisa < 0x800000)
		{
			mantisa <<= 1;
			exponent--;
		}

	}
	
	mantisa |= 0x800000;

ret:
	return {mantisa, exponent, sign};
}

constexpr uint32_t pack(const Content& c)
{
	uint32_t sign     = c.sign;
	uint32_t exponent = c.exponent;
	uint32_t mantisa  = c.mantisa;

	//arbitrary border for overflow detection
	if(exponent > 4096)
		return (sign << 31) | 0;
	if(c.exponent == 0xFF + 1024 - 127)
	{
		return QNAN;
	}
	if(exponent > 1024 + 127)
	{
		return (sign << 31) | (INF - 1);
	}

	if(exponent <= 1024 - 127)
		exponent--;

	//convert from extended precision normal to subnormal
	while(exponent < 1024 - 127
	   && mantisa != 0)
	{
		exponent++;
		mantisa >>= 1;
	}

	//if no implicit one, return 0 because extended precision 
	if(mantisa == 0) 
		return (sign << 31);
	
	return (( sign                     & 0x1)  << 31)
	     | (((exponent - (1024 - 127)) & 0xFF) << 23)
	     | (( mantisa & 0x7FFFFF    ));

}

/*
3.4 binary interchange format encodings

The representation r of the floating-point datum, and value v of the floating-point datum represented,
are inferred from the constituent fields as follows:
a) if EXPONENT = 2^w - 1 (all ones) and MANTISA != 0 then r is QNAN and SNAN regardless of SIGN
b) if EXPONENT = 2^w - 1 and MANTISA = 0 then (if SIGN = 0 r is infinity else r is negative infitiy)
[...]
e) if EXPONENT = 0 and MANTISA = 0 then (if SIGN = 0 r is +0 else r is -0)
*/
constexpr bool is_NaN(const uint32_t a)
{
	return ((a >> 23) & 0xFF) == 0xFF
	       && ((a & 0x7FFFFF) != 0); 
}
constexpr bool is_sNaN(const uint32_t a)
{
	return ((a >> 23) & 0xFF) == 0xFF
	       && ((a & 0x7FFFFF) != 0)
	       && ((a & 0x400000) == 0);
}
constexpr bool is_Fin(const uint32_t a)
{
	return ((a >> 23) & 0xFF) != 0xFF;
}
constexpr bool is_Inf(const uint32_t a)
{
	return ((a >> 23) & 0xFF) == 0xFF
	       && ((a & 0x7FFFFF) == 0); 
}
constexpr bool is_Zero(const uint32_t a)
{
	return (a & 0x7FFFFFFF) == 0;
}

constexpr bool is_Neg(const uint32_t a)
{
	return (a & 0x80000000);
}


/*
If x is the negative number of least magnitude in x’s format, nextUp(x) is −0. 
 nextUp(+-0) is the positive number of least magnitude in x’s format. 
 nextUp(+inf) is +inf
 nextUp(−inf) is the finite negative number largest in magnitude.
 When x is NaN, then the result is according to 6.2. 
 nextUp(x) is quiet except for sNaNs.
*/

constexpr uint32_t next_Up(const uint32_t a)
{
	if(is_NaN(a))
		return a;
	if(is_Zero(a))
		return 0x1;

	if(is_Neg(a))
		return a - 1;
	else
	{
		if(is_Inf(a))
			return a;
	
		return a + 1;
	}
}

/*
	as defined by the standard
*/
constexpr uint32_t next_Down(const uint32_t a)
{
	return neg(next_Up(neg(a)));
}

/*
	for justification of comparison functions, see standard
*/
constexpr bool qeq(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a) || is_NaN(b))
		return false;

	if(is_Zero(a) && is_Zero(b)) return true;

	return a == b;
}

constexpr bool qne(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a) || is_NaN(b))
		return false;
	
	if(is_Zero(a) && is_Zero(b)) return false;

	return a != b;
}

constexpr bool qgr(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a) || is_NaN(b))
		return false;

	const bool sign_a = a >> 31;
	const bool sign_b = b >> 31;

	//both positive
	if(!sign_a && !sign_b)
		return a > b;
	//both negative
	if(sign_a && sign_b)
		return b > a;

	if(sign_a && !sign_b)
		return true;
	else
	//if(sign_b && !sign_a) not needed but left as comment for readability
		return false;
}


constexpr bool qge(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a) || is_NaN(b))
		return false;

	const bool sign_a = a >> 31;
	const bool sign_b = b >> 31;

	//both positive
	if(!sign_a && !sign_b)
		return a >= b;
	//both negative
	if(sign_a && sign_b)
		return b >= a;

	if(sign_a && !sign_b)
		return true;
	else
	//if(sign_b && !sign_a) not needed but left as comment for readability
		return false;
}

constexpr bool qls(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a) || is_NaN(b))
		return false;

	const bool sign_a = a >> 31;
	const bool sign_b = b >> 31;

	//both positive
	if(!sign_a && !sign_b)
		return a < b;
	//both negative
	if(sign_a && sign_b)
		return b < a;

	if(sign_a && !sign_b)
		return false;
	else
	//if(sign_b && !sign_a) not needed but left as comment for readability
		return true;
}

constexpr bool qle(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a) || is_NaN(b))
		return false;

	const bool sign_a = a >> 31;
	const bool sign_b = b >> 31;

	//both positive
	if(!sign_a && !sign_b)
		return a <= b;
	//both negative
	if(sign_a && sign_b)
		return b <= a;

	if(sign_a && !sign_b)
		return false;
	else
	//if(sign_b && !sign_a) not needed but left as comment for readability
		return true;
}

constexpr uint32_t mul(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a))
		return a;
	if(is_NaN(b))
		return b;
	
	if((is_Inf(a) && is_Zero(b))
	|| (is_Zero(a) && is_Inf(b)))
		return QNAN;
		
	//sign remains only when the sign is different
	const uint8_t sign = (a ^ b) >> 31; 
	if(is_Inf(a)
	|| is_Inf(b))
	{
		return static_cast<uint32_t>(sign) << 31 | INF;
	}
	
	const Content f = unpack(a);
	const Content s = unpack(b);

	
	/*
		first exponents are debiasaed, then added, then sum is biased again
		mantisas are multiplied like integers and lower bits are discarded
	*/
	uint16_t exponent =   f.exponent
	                    + s.exponent
	                    - BIAS;
	uint64_t mantisa =   static_cast<uint64_t>(f.mantisa) 
	                   * static_cast<uint64_t>(s.mantisa);
	mantisa >>= 23;

	//infinite loop prevention
	if(mantisa == 0)
		goto mul_ret;

	//while mantisa not in range, adjust
	while(mantisa > 0xFFFFFF) //24 bits
	{
		exponent++;
		mantisa >>= 1;
	} 
	while(mantisa < 0x800000) //24 bits
	{
		exponent--;
		mantisa <<= 1;
	} 

mul_ret:
	const uint32_t t = pack({static_cast<uint32_t>(mantisa), exponent, sign});

	//preveention of wrong result due to representation, same later
	if(is_NaN(t))
		return (static_cast<uint32_t>(sign) << 31) | (INF - 1);

	return t;
}

constexpr uint32_t div(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a))
		return a;
	if(is_NaN(b))
		return b;
	
	const uint32_t sign = (a ^ b) & 0x80000000;
	
	IF_A_B_BOTH(is_Inf(a), is_Inf(b), 
	            sign | INF, sign | 0, QNAN);
	IF_A_B_BOTH(is_Zero(a), is_Zero(b), 
	            sign | 0, sign | INF, QNAN);

	Content f = unpack(a);
	Content s = unpack(b);

	f.sign ^= s.sign;

	/*
		shifted 40 to increase precision of result
		later res is shifted 40 to the right
	*/
	uint64_t num = static_cast<uint64_t>(f.mantisa) << 40;
	uint64_t div = s.mantisa;
	uint64_t res = num / div;

	res <<= 23;

	/*
		subtracting exponents so BIAS has to be added
	*/
	f.exponent = f.exponent - s.exponent + BIAS;

	if(res == 0)
		return f.sign | (INF - 1);

	while((res & 0x8000000000000000) == 0)
	{
		res <<= 1;
		f.exponent--;
	}
	
	f.mantisa = res >> 40;
	const uint32_t t = pack(f);
	if(is_NaN(t))
		return (static_cast<uint32_t>(f.sign) << 31) | (INF - 1);

	return t;
}

constexpr uint32_t add(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a))
		return a;
	if(is_NaN(b))
		return b;

	if(is_Neg(b))
		return sub(a, neg(b));
	if(is_Neg(a))
		return sub(b, neg(a));

	//they are positive
	if(is_Inf(a)
	|| is_Inf(b))
	{
		return INF - 1;
	}

	Content f = unpack(a);
	Content s = unpack(b);


	const uint8_t sign = 0;

	//adjust smaller number and simply add their mantisa
	//then adjust result
	uint16_t exponent;
	uint32_t to_shift;
	uint32_t mantisa;
	if(f.exponent < s.exponent)
	{
		exponent = s.exponent; 
		to_shift = s.exponent - f.exponent;
		mantisa  = (f.mantisa >> to_shift) + s.mantisa;
		if(to_shift >= 32)
			mantisa = s.mantisa;
	}
	else
	{
		exponent = f.exponent; 
		to_shift = f.exponent - s.exponent;
		mantisa  = (s.mantisa >> to_shift) + f.mantisa;
		if(to_shift >= 32)
			mantisa = f.mantisa;
	}
	while(mantisa > 0xFFFFFF) //24 bits
	{
		exponent++;
		mantisa >>= 1;
	} 

	const uint32_t t = pack({mantisa, exponent, sign});
	if(is_NaN(t))
		return (static_cast<uint32_t>(sign) << 31) | (INF - 1);

	return t;
}

//baically same logic as add
constexpr uint32_t sub(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a))
		return a;
	if(is_NaN(b))
		return b;

	if(is_Neg(b))
		return add(a, neg(b));
	if(is_Neg(a))
		return neg(add(neg(a), b));


	IF_A_B_BOTH(is_Zero(a), is_Zero(b),
	            neg(b), a, 0);
	IF_A_B_BOTH(is_Inf(a), is_Inf(b),
	            INF, neg(INF), QNAN);


	Content f = unpack(a);
	Content s = unpack(b);

	const uint8_t sign = qgr(b, a);

	uint16_t exponent; 
	uint32_t to_shift; 
	uint64_t fmantisa = f.mantisa; 
	uint64_t smantisa = s.mantisa; 
	uint64_t mantisa;

	fmantisa <<= 40;
	smantisa <<= 40;

	if(!sign)
	{
		exponent = f.exponent;
		to_shift = f.exponent - s.exponent;
		mantisa  = fmantisa - (smantisa >> to_shift);

		if(to_shift >= 64)
			mantisa = fmantisa - 1;
		
	}
	else
	{
		exponent = s.exponent;
		to_shift = s.exponent - f.exponent;
		mantisa  = smantisa - (fmantisa >> to_shift);
		if(to_shift >= 64)
			mantisa = smantisa - 1;
	}

	if(mantisa == 0)
		return 0;

	while(mantisa > (0xFFFFFFUL << 40)) //24 bits
	{
		exponent++;
		mantisa >>= 1;
	} 
	while(mantisa < (0x800000UL << 40)) //24 bits
	{
		exponent--;
		mantisa <<= 1;
	} 
	
	if(mantisa == 0)
		return 0;
	mantisa >>= 40;
	
	const uint32_t t = pack({static_cast<uint32_t>(mantisa), exponent, sign});
	if(is_NaN(t))
		return (static_cast<uint32_t>(sign) << 31) | (INF - 1);

	return t;
}

constexpr uint32_t neg(const uint32_t a)
{
	if(is_NaN(a))
		return a;
	
	return a ^ (0x80000000);
}

constexpr uint32_t abs(const uint32_t a)
{
	if(is_NaN(a))
		return a;
	
	return a & ~(0x80000000);
}

constexpr uint32_t whl(const uint32_t a)
{
	if(is_Inf(a)
	|| is_NaN(a)
	|| is_Zero(a))
		return a;

	Content cont = unpack(a);

	const int32_t to_shift = cont.exponent - BIAS;
	if(to_shift < 0)
		return 0;
	if(to_shift > 23)
		return a;
	cont.mantisa >>= 23 - (cont.exponent - BIAS);
	cont.mantisa <<= 23 - (cont.exponent - BIAS);

	return pack(cont);
}

constexpr uint32_t min_Val(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a))
		return b;
	if(is_NaN(b))
		return a;

	if(qgr(a, b))
		return b;
	else
		return a;
}
constexpr uint32_t max_Val(const uint32_t a, const uint32_t b)
{
	if(is_NaN(a))
		return b;
	if(is_NaN(b))
		return a;
	
	if(qgr(a, b))
		return a;
	else
		return b;
}

constexpr uint32_t scl(const uint32_t a, const int32_t N)
{
	if(is_Zero(a))
		return a;
	if(is_Inf(a))
		return a;
	if(is_NaN(a))
		return a;

	if(N >  1023)
		return INF - 1;
	if(N < -1023)
		return neg(INF) - 1;


	Content f = unpack(a);
	f.exponent += N;
	const uint32_t t = pack(f);
	return is_Inf(t) ? t - 1 : t;
}

constexpr int32_t log(const uint32_t a)
{
	//check for NaN not needed, since NaN is negative
	if(is_Zero(a)
	|| is_Inf(a))
		return INT_MIN; //value outside of range 

	Content cont = unpack(a);
	return cont.exponent - BIAS;
}

constexpr int32_t fl_To_Int(const uint32_t a)
{
//returns INT_MIN to be compatible with c++

	if(is_Inf(a))
	{
		//is negative
//		if(fl & 0x80000000)
			return INT_MIN;
		
		return INT_MAX;
	}

	Content cont = unpack(a);

//	std::cout << '\n';
//	std::cout << BIAS + 32 << '\n';
//	std::cout << +cont.sign << ' ' << ' ' << cont.exponent << ' ' << cont.mantisa << '\n'; 
	
	if(cont.exponent < BIAS)
		cont.mantisa = 0;
	else if(cont.exponent >= BIAS + 31)
	{
		//is negative
//		if(fl & 0x80000000)
			return INT_MIN;
		
		return INT_MAX;
	}
	else if(cont.exponent > BIAS + 23)
	{
		cont.mantisa <<= cont.exponent - BIAS - 23;
	}
	else
		cont.mantisa >>= 23 - (cont.exponent - BIAS);

	return cont.sign ? -cont.mantisa : cont.mantisa;
}

constexpr uint32_t cps(const uint32_t a, const uint32_t b)
{
	uint32_t ret  = a & 0x7FFFFFFF;
	         ret |= b & 0x80000000;	
	return ret;
}

constexpr uint32_t int_To_Fl(int32_t i)
{
	if(i == 0) return 0;

	uint32_t sign = 0;
	if(i < 0)
	{
		sign = 0x80000000;
		i = -i;
	}

	uint32_t mantisa  = i;
	uint32_t exponent = 127 + 23;

	while(mantisa > 0xFFFFFF)
	{
		mantisa >>= 1;
		exponent++;
	}
	while(mantisa < 0x800000)
	{
		mantisa <<= 1;
		exponent--;
	}

	exponent <<= 23;
	mantisa &= 0x7FFFFF;

	return sign | exponent | mantisa;
}

constexpr bool is_Sbn(const uint32_t a)
{
	return ((a >> 23) & 0xFF) == 0;
}

constexpr bool is_Nrm(const uint32_t a)
{
	if(is_NaN(a)
	|| is_Inf(a))
		return false;

	return !is_Sbn(a);
}

constexpr uint32_t sqt(const uint32_t a)
{
	if(is_Zero(a)) return a;
	if(is_Neg(a)) return QNAN;
	if(is_Inf(a)) return INF;
	

	Content f = unpack(a);

	/*
	sqrt(2^exp * mant)
	sqrt(2^exp) * sqrt(mant) 
	2^(exp / 2) * sqrt(mant)

	sqrt(mant) can be computed as integer

	in case of odd number we have to multiply mantisa by 2 first, then take square root
	because 2^(odd / 2) will lead to 2^(even + 1/2) = 2^even * sqrt(2) 
	sqrt(2 * mantisa)  = sqrt(mantisa) * sqrt(2)
	*/

	f.exponent -= BIAS;
	const bool odd = f.exponent & 0b1;
	f.exponent  = static_cast<int16_t>(f.exponent) >> 1;

	f.exponent += BIAS;

	uint64_t tmp = f.mantisa;
	tmp = std::sqrt(tmp << (39 + odd));
	tmp >>= 8;
	f.mantisa = tmp;

	return pack(f);

}

constexpr cls float_Class(const uint32_t a)
{
	using enum cls;
	if(is_sNaN(a))
		return snan;
	if(is_NaN(a))
		return qnan;
	if(is_Inf(a))
	{
		if(is_Neg(a))
			return ninf;
		else
			return pinf;
	}
	if(is_Zero(a))
	{
		if(is_Neg(a))
			return nzro;
		else
			return pzro;
	}
	if(is_Nrm(a))
	{
		if(is_Neg(a))
			return nnrm;
		else
			return pnrm;
	}
	if(is_Neg(a))
		return nsbn;
	else
		return psbn;
	
};

/*
	section 5.10 details of totalOrder predicate

	i am not sure if i correctly understand what should happen when standard doesnt specify what should be the return value
	eg is_NaN(a) && !is_NaN(b) && !is_Neg(a)
	i assume it should return true
*/
constexpr bool total_Order(const uint32_t a, const uint32_t b)
{
	//d)
	if(is_NaN(a))
	{
		//3)
		if(is_NaN(b))
		{
			//i)
			if(is_Neg(a))
				return true;

			//ii)
			//???

			//iii)
			//?????

			return false;
		}
		//1)
		else
		{
			//1)
			if(is_Neg(a))
				return true;
	
			return false;
		}
	}
	if(is_NaN(b))
	{
		//2)
		if(is_Neg(b))
			return false;

		return true;
	}


	//a)
	if(qls(a, b))
		return true;
	//b)
	if(qgr(a, b))
		return false;

	//c)
	if(a == b)
		return true;
	//1) and 2)
	//if(is_Zero(a) && is_Zero(b)) ; not needed since if qeq() == true but a != b then they must be zeros of different sign
	{
		if(is_Neg(a))
			return true;

		return false;
	}

}


int main()
{
	std::cout << std::hex;
	std::fesetround(FE_TOWARDZERO);
	
	union{
		float    f;
		uint32_t i;
	} n1;
	union{
		float    f;
		uint32_t i;
	} n2;

	n1.i = 0x80000000;
	n2.i = 0x800FFFFF;
	for(n1.i = 0x00000000; n1.i < 0x7F7FFFFF; n1.i += 0xFFFFF)
	for(n2.i = 0x00000000; n2.i < 0x7F7FFFFF; n2.i += 0xFFFFF)
//	for(n1.i = 0x80000000; n1.i < 0xFF7FFFFF; n1.i += 0xFFFFF)
//	for(n2.i = 0x80000000; n2.i < 0xFF7FFFFF; n2.i += 0xFFFFF)
{
	std::cout << n1.i << ' ' << n2.i << '\n';
	std::cout << n1.f << ' ' << n2.f << '\n';
	
	float f = n1.f;
	float s = n2.f;

	Content conf = unpack(*reinterpret_cast<uint32_t*>(&f));
	Content cons = unpack(*reinterpret_cast<uint32_t*>(&s));

	uint32_t valf = pack(conf);
	uint32_t vals = pack(cons);
	uint32_t res;

#define GNRL 1
#define LOG2 2
#define ARIT 3
#define SIGN 4
#define COMP 5
#define NCVR 6
#define NCGN 7
#define NCFL 8
#define OP NCGN 

	std::cout << "c  ";
	CHECK(valf, f);
	CHECK(vals, s);

#if OP == GNRL	
	std::cout << "\nnp ";
	CHECK(next_Up(valf), std::nextafter(f, INFINITY));
	CHECK(next_Up(vals), std::nextafter(s, INFINITY));
	
	std::cout << "\nnd ";
	CHECK(next_Down(valf), std::nextafter(f, -INFINITY));
	CHECK(next_Down(vals), std::nextafter(s, -INFINITY));

	std::cout << "\nmnv";
	CHECK(min_Val(valf, vals), std::fmin(f, s));
	CHECK(min_Val(vals, valf), std::fmin(s, f));
	std::cout << "\nmxv";
	CHECK(max_Val(valf, vals), std::fmax(f, s));
	CHECK(max_Val(vals, valf), std::fmax(s, f));

#elif OP == LOG2
	std::cout << "\ns2 ";
	CHECK(scl(valf, fl_To_Int(vals)), std::scalbn(f, static_cast<int>(s)));
	CHECK(scl(vals, fl_To_Int(valf)), std::scalbn(s, static_cast<int>(f)));

	std::cout << "\nl2 ";
	CHECK_INT(log(valf), static_cast<int>(std::logb(f)));
	CHECK_INT(log(vals), static_cast<int>(std::logb(s)));

#elif OP == SIGN
	std::cout << "\nneg";
	CHECK(neg(valf), -f);
	CHECK(neg(vals), -s);

	std::cout << "\nabs";
	CHECK(abs(valf), std::abs(f));
	CHECK(abs(vals), std::abs(s));

	std::cout << "\ncps";
	CHECK(cps(valf, vals), std::copysign(f, s));
	CHECK(cps(vals, valf), std::copysign(s, f));

#elif OP == COMP

	std::cout << "\nqeq";
	CHECK_INT(qeq(valf, vals), f == s);
	CHECK_INT(qeq(vals, valf), s == f);
	
	std::cout << "\nqne";
	CHECK_INT(qne(valf, vals), f != s);
	CHECK_INT(qne(vals, valf), s != f);
	
	std::cout << "\nqgr";
	CHECK_INT(qgr(valf, vals), f >  s);
	CHECK_INT(qgr(vals, valf), s >  f);
	
	std::cout << "\nqge";
	CHECK_INT(qge(valf, vals), f >= s);
	CHECK_INT(qge(vals, valf), s >= f);
	
	std::cout << "\nqls";
	CHECK_INT(qls(valf, vals), f <  s);
	CHECK_INT(qls(vals, valf), s <  f);
	
	std::cout << "\nqle";
	CHECK_INT(qle(valf, vals), f <= s);
	CHECK_INT(qle(vals, valf), s <= f);



#elif OP == ARIT
	std::cout << "\nmul";
	CHECK(mul(valf, vals), f * s);
	CHECK(mul(vals, valf), s * f);
	
	std::cout << "\ndiv";
	CHECK(div(valf, vals), f / s);
	CHECK(div(vals, valf), s / f);

	std::cout << "\nadd";
	CHECK(add(valf, vals), f + s);
	CHECK(add(vals, valf), s + f);
	
	std::cout << "\nsub";
	CHECK(sub(valf, vals), f - s);
	CHECK(sub(vals, valf), s - f);
	
	std::cout << "\nsqt";
	CHECK(sqt(valf), std::sqrt(f));
	CHECK(sqt(vals), std::sqrt(s));
	
	std::cout << "\ntoi";
	CHECK_INT(fl_To_Int(valf), static_cast<int32_t>(f));
	CHECK_INT(fl_To_Int(vals), static_cast<int32_t>(s));
	
	std::cout << "\nfri";
	CHECK(int_To_Fl(fl_To_Int(valf)), static_cast<float>(static_cast<int32_t>(f)));
	CHECK(int_To_Fl(fl_To_Int(vals)), static_cast<float>(static_cast<int32_t>(s)));


#elif OP == NCVR

	std::cout << "\nIEEE 1985\n";
	std::cout << is_754_Version_1985() << '\n';
	std::cout << "IEEE 2008\n";
	std::cout << is_754_Version_2008() << '\n';

#elif OP == NCGN

	std::cout << "\nneg";
	CHECK_INT(is_Neg(valf), (f < 0));
	CHECK_INT(is_Neg(vals), (s < 0));

	std::cout << "\nnrm";
	CHECK_INT(is_Nrm(valf), std::isnormal(f));
	CHECK_INT(is_Nrm(vals), std::isnormal(s));
	
	std::cout << "\nsbn";
	CHECK_INT(is_Sbn(valf), !std::isnormal(f));
	CHECK_INT(is_Sbn(vals), !std::isnormal(s));

	std::cout << "\ninf";
	CHECK_INT(is_Inf(valf), std::isinf(f));
	CHECK_INT(is_Inf(vals), std::isinf(s));

	std::cout << "\nnan";
	CHECK_INT(is_NaN(valf), std::isnan(f));
	CHECK_INT(is_NaN(vals), std::isnan(s));
#endif 

	std::cout << "\n\n";
}
}

