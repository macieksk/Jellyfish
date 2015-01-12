/*
 * seedmod.hpp
 *
 *  Created on: Dec 16, 2014
 *      Author: Maciek Sykulski <macieksk@gmail.com>
 */

#ifndef SEEDMOD_HPP_
#define SEEDMOD_HPP_

#include <string.h>

#include <iterator>

namespace jellyfish {
template <typename base_type,
		  class MerTypeOut>
class mer_shift_left_output_iterator:
		public std::iterator<std::output_iterator_tag,MerTypeOut>
{
private:
	MerTypeOut & mer_;

public:
	typedef mer_shift_left_output_iterator<base_type,MerTypeOut> this_type;

	mer_shift_left_output_iterator(MerTypeOut & mer)
			:mer_(mer){}

	this_type& operator=(const base_type & c) {
		mer_.shift_left((const int &)c);
		return *this;
	}

	this_type& operator*() {return *this;}
	this_type& operator++(){return *this;} //prefix increment
};


template <typename base_type,
		  class MerTypeOut>
class mer_explicit_shift_left_output_iterator:
		public std::iterator<std::output_iterator_tag,MerTypeOut>
{
private:
	MerTypeOut & mer_;

public:
	typedef mer_explicit_shift_left_output_iterator<base_type,MerTypeOut> this_type;

	mer_explicit_shift_left_output_iterator(MerTypeOut & mer)
			:mer_(mer){}

	this_type& operator=(const base_type & c) {
		mer_.shift_left((const int &)c);
		return *this;
	}

	this_type& operator*() {return *this;}
	this_type& operator++(){return *this;} //prefix increment
};


}

namespace seedmod {
  	//This is an OutputIterator
    //operator= <--> ret_m.shift_left(*c);
	template <typename InputType,
			  class OutMerIterator,
		      char MATCH_CHR = '#',
		      char DEL_CHR = '^',
			  char INS_CHR = 'v'>
	class SpacedSeedForIndexSquasherIterator
	{
	private:
		const char * seed_;
		OutMerIterator & oit_;

	public:
		typedef SpacedSeedForIndexSquasherIterator<InputType,OutMerIterator,
				MATCH_CHR,DEL_CHR,INS_CHR> this_type;

		SpacedSeedForIndexSquasherIterator(const char * seed, OutMerIterator & oit)
			:seed_(seed),oit_(oit){}

		this_type& operator=(const InputType c) {
			do{
				if(DEL_CHR==*seed_){
					++seed_;
					continue;
				}else if(MATCH_CHR==*seed_){
					*oit_ = c;
				}else if(INS_CHR==*seed_){
				}
				//space
				break;
			}while (true);
			++seed_;
			return *this;
		}

		this_type& operator*() {return *this;}
		this_type& operator++(){return *this;} //prefix increment

		bool at_end() const{
			return seed_==0;
		}

	public:
		static size_t weight(const char * s){
			size_t w=0;
			for (;*s;++s)
				if (*s==MATCH_CHR) ++w;
			return w;
		}
		static size_t span(const char * s){
			size_t w=0;
			for (;*s;++s)
				if (*s!=INS_CHR) ++w;
			return w;
		}
		static bool isfullseed(const char * s){
					for (;*s;++s)
						if (*s!=MATCH_CHR)
							return false;
					return true;
		}
	};

	//ForReadSquasher is the same as ForIndexSquasher
	//but with DEL_CHR, INS_CHR swapped
	template <typename InputType,
				class OutMerIterator,
			      char MATCH_CHR = '#',
			      char DEL_CHR = '^',
				  char INS_CHR = 'v'>
	class SpacedSeedForReadSquasherIterator:
			public SpacedSeedForIndexSquasherIterator<InputType,OutMerIterator,
													MATCH_CHR,INS_CHR,DEL_CHR>
	{
	public:
		typedef SpacedSeedForReadSquasherIterator<InputType,OutMerIterator,
				MATCH_CHR,DEL_CHR,INS_CHR> this_type;
		typedef SpacedSeedForIndexSquasherIterator<InputType,OutMerIterator,
				MATCH_CHR,INS_CHR,DEL_CHR> super;

		SpacedSeedForReadSquasherIterator(const char * seed, OutMerIterator & oit)
			:super(seed,oit){}
	};

}


namespace kraken {

template<typename base_type, typename OutputIterator>
static OutputIterator to_codes(const size_t k, const base_type & w, OutputIterator it) {
	static const base_type c3 = (base_type)0x3;
	int shift  = (k<<1) - 2; // Number of bits to shift to get base
  for( ; shift >= 0; shift -= 2, ++it)
      *it = (w >> shift) & c3;
  return it;
}



template<typename base_type>
class kmer_shift_left_output_iterator:
		public std::iterator<std::output_iterator_tag,base_type>
{
private:
	base_type & kmer_;
	static const base_type c3 = (base_type)0x3;

public:
	typedef kmer_shift_left_output_iterator<base_type> this_type;

	kmer_shift_left_output_iterator(base_type & kmer)
			:kmer_(kmer){}

	this_type& operator=(const base_type & c) {
		//mer_.shift_left(c);
		kmer_<<=2;
		kmer_ |= c & c3;
		return *this;
	}

	this_type& operator*() {return *this;}
	this_type& operator++(){return *this;} //prefix increment
};

//Chaining input/output iterators to squash seed
template<typename base_type>
static void squash_kmer_for_read(const char * seed, size_t seed_len,
								 const base_type & fmer, base_type & ret_m){
    typedef kraken::kmer_shift_left_output_iterator<base_type>
    									mer_sleft_oiter;
    typedef seedmod::SpacedSeedForReadSquasherIterator<base_type,
      		  	  	  	  	  	  	  	  	  mer_sleft_oiter> seed_read_squasher_iter_type;
  	mer_sleft_oiter meroiter(ret_m);
  	seed_read_squasher_iter_type seed_squash_it(seed,meroiter);
  	kraken::to_codes(seed_len,fmer,seed_squash_it);
}

template<typename base_type>
static void squash_kmer_for_index(const char * seed, size_t seed_len,
								  const base_type & fmer, base_type & ret_m){
    typedef kraken::kmer_shift_left_output_iterator<base_type>
    									mer_sleft_oiter;
	typedef seedmod::SpacedSeedForIndexSquasherIterator<base_type,
          		  	  	  	  	  	  	  	  	  mer_sleft_oiter> seed_index_squasher_iter_type;
    mer_sleft_oiter meroiter(ret_m);
    seed_index_squasher_iter_type seed_squash_it(seed,meroiter);
    kraken::to_codes(seed_len,fmer,seed_squash_it);
}

}






#endif /* SEEDMOD_HPP_ */


//   typedef SpacedSeed<'#'> spaced_seed;
//
//#include <boost/iterator/iterator_facade.hpp>
//
//   namespace impl {
//	   template <class Value>
//	   class node_iter
//		 : public boost::iterator_facade<
//			   node_iter<Value>
//			 , Value
//			 , boost::forward_traversal_tag
//		   >
//	   {
//		public:
//		   node_iter()
//			 : m_node(0) {}
//
//		   explicit node_iter(Value* p)
//			 : m_node(p) {}
//
//		   template <class OtherValue>
//		   node_iter(node_iter<OtherValue> const& other)
//			 : m_node(other.m_node) {}
//
//		private:
//		   friend class boost::iterator_core_access;
//		   template <class> friend class node_iter;
//
//		   template <class OtherValue>
//		   bool equal(node_iter<OtherValue> const& other) const
//		   {
//			   return this->m_node == other.m_node;
//		   }
//
//		   void increment()
//		   { m_node = m_node->next(); }
//
//		   Value& dereference() const
//		   { return *m_node; }
//
//		   Value* m_node;
//	   };
//   }
//   typedef impl::node_iter<spaced_seed> node_iterator;
//   typedef impl::node_iter<spaced_seed const> node_const_iterator;

