#ifndef __MEM_OPER_DEFS_H__
#define __MEM_OPER_DEFS_H__


#define MEM_UNARY_OPER(CLASS, OP)      	       	       	       	       	    \
template<typename T>							    \
void CLASS<T>::OP(T* a_o, const T* a_i, size_t n, StreamT stream){	    \
    Comp_unary< MATH_##OP>(a_o, a_i, n, stream);			    \
}

#define MEM_UNARY_MASKED_OPER(CLASS, OP)       	       	       	       	    \
template<typename T>							    \
void CLASS<T>::								    \
OP(T* a_o, const T* a_i, const T* a_mask,				    \
   size_t n, StreamT stream)						    \
{									    \
    Comp_unary< MATH_##OP>(a_o, a_i, a_mask, n, stream);		    \
}

#define MEM_UNARY_OPER_I(CLASS, OP)					    \
template<typename T>							    \
void CLASS<T>::OP##_I(T* a_o, size_t n, StreamT stream){		    \
    Comp_unary_I<MATH_##OP>(a_o, n, stream);				    \
}

#define MEM_UNARY_MASKED_OPER_I(CLASS, OP)				    \
template<typename T>							    \
void CLASS<T>::								    \
OP##_I(T* a_o, const T* a_mask, size_t n, StreamT stream)		    \
{									    \
    Comp_unary_I<MATH_##OP>(a_o, a_mask, n, stream);                        \
}

#define MEM_UNARY_OPERS(CLASS, OP)					    \
    MEM_UNARY_OPER(CLASS, OP)						    \
    MEM_UNARY_MASKED_OPER(CLASS, OP)					    \
    MEM_UNARY_OPER_I(CLASS, OP)                                             \
    MEM_UNARY_MASKED_OPER_I(CLASS, OP)

#define MEM_BINARY_OPER_NOC(CLASS, OP) 	       	       	       	       	    \
template<typename T>							    \
void CLASS<T>::OP(T* a_o, const T* a_i, const T& c,			    \
		  size_t n, StreamT stream, bool onDev)			    \
{									    \
    binaryC<MATH_##OP>(a_o, a_i, c, n, stream, onDev);			    \
}

#define MEM_BINARY_MASKED_OPER_NOC(CLASS, OP)	       	       	       	    \
template<typename T>							    \
void CLASS<T>::OP(T* a_o, const T* a_i, const T& c, const T* a_mask,	    \
		  size_t n, StreamT stream, bool onDev)			    \
{									    \
    binaryC<MATH_##OP>(a_o, a_i, c, a_mask, n, stream, onDev);		    \
}

#define MEM_BINARY_OPER_NOC_I(CLASS, OP)                                    \
template<typename T>							    \
void CLASS<T>::OP##_I(T* a_o, const T& c,				    \
		      size_t n, StreamT stream, bool onDev)		    \
{									    \
    binaryC_I<MATH_##OP>(a_o, c, n, stream, onDev);			    \
}

#define MEM_BINARY_MASKED_OPER_NOC_I(CLASS, OP)				    \
template<typename T>							    \
void CLASS<T>::OP##_I(T* a_o, const T& c, const T* a_mask,		    \
		      size_t n, StreamT stream, bool onDev)		    \
{									    \
    binaryC_I<MATH_##OP>(a_o, c, a_mask, n, stream, onDev);		    \
}

#define MEM_BINARY_OPERC(CLASS, OP)					    \
template<typename T>							    \
void CLASS<T>::OP##C(T* a_o, const T* a_i, const T& c,			    \
		      size_t n, StreamT stream, bool onDev)		    \
{									    \
    binaryC<MATH_##OP>(a_o, a_i, c, n, stream, onDev);			    \
}

#define MEM_BINARY_MASKED_OPERC(CLASS, OP)				    \
template<typename T>							    \
void CLASS<T>::OP##C(T* a_o, const T* a_i, const T& c, const T* a_mask,	    \
		      size_t n, StreamT stream, bool onDev)		    \
{									    \
    binaryC<MATH_##OP>(a_o, a_i, c, a_mask, n, stream, onDev);		    \
}

#define MEM_BINARY_OPERC_I(CLASS, OP)					    \
template<typename T>							    \
void CLASS<T>::OP##C_I(T* a_o, const T& c,				    \
			  size_t n, StreamT stream, bool onDev)		    \
{									    \
    binaryC_I<MATH_##OP>(a_o, c, n, stream, onDev);			    \
}

#define MEM_BINARY_MASKED_OPERC_I(CLASS, OP)				    \
template<typename T>							    \
void CLASS<T>::OP##C_I(T* a_o, const T& c, const T* a_mask,		    \
			  size_t n, StreamT stream, bool onDev)		    \
{									    \
    binaryC_I<MATH_##OP>(a_o, c, a_mask, n, stream, onDev);		    \
}

#define MEM_BINARY_OPER(CLASS, OP)					    \
template<typename T>							    \
void CLASS<T>::OP(T* a_o, const T* a_i, const T* a_i1,			    \
		      size_t n, StreamT stream)				    \
{									    \
    binary<MATH_##OP>(a_o, a_i, a_i1, n, stream);			    \
}

#define MEM_BINARY_MASKED_OPER(CLASS, OP)				    \
template<typename T>							    \
void CLASS<T>::OP(T* a_o, const T* a_i, const T* a_i1, const T* a_mask,	    \
		      size_t n, StreamT stream)				    \
{									    \
    binary<MATH_##OP>(a_o, a_i, a_i1, a_mask, n, stream);		    \
}

#define MEM_BINARY_OPER_I(CLASS, OP)					    \
template<typename T>							    \
void CLASS<T>::OP##_I(T* a_o, const T* a_i,				    \
			  size_t n, StreamT stream)			    \
{									    \
    binary_I<MATH_##OP>(a_o, a_i, n, stream);				    \
}

#define MEM_BINARY_MASKED_OPER_I(CLASS, OP)				    \
template<typename T>							    \
void CLASS<T>::OP##_I(T* a_o, const T* a_i, const T* a_mask,		    \
			  size_t n, StreamT stream)			    \
{									    \
    binary_I<MATH_##OP>(a_o, a_i, a_mask, n, stream);			    \
}

#define MEM_BINARY_OPERS(CLASS, OP)					    \
MEM_BINARY_OPERC(CLASS, OP)						    \
MEM_BINARY_MASKED_OPERC(CLASS, OP)					    \
MEM_BINARY_OPERC_I(CLASS, OP)						    \
MEM_BINARY_MASKED_OPERC_I(CLASS, OP)					    \
MEM_BINARY_OPER_NOC(CLASS, OP)						    \
MEM_BINARY_MASKED_OPER_NOC(CLASS, OP)					    \
MEM_BINARY_OPER_NOC_I(CLASS, OP)					    \
MEM_BINARY_MASKED_OPER_NOC_I(CLASS, OP)					    \
MEM_BINARY_OPER(CLASS, OP)						    \
MEM_BINARY_MASKED_OPER(CLASS, OP)					    \
MEM_BINARY_OPER_I(CLASS, OP)                                                \
MEM_BINARY_MASKED_OPER_I(CLASS, OP)

#endif // __MEM_OPER_DEFS_H__



