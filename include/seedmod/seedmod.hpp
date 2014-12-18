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
template <class MerTypeOut>
class mer_shift_left_output_iterator:
		public std::iterator<std::output_iterator_tag,MerTypeOut>
{
private:
	MerTypeOut & mer_;

public:
	typedef mer_shift_left_output_iterator<MerTypeOut> this_type;

	mer_shift_left_output_iterator(MerTypeOut & mer)
			:mer_(mer){}

	this_type& operator=(const char c) {
		mer_.shift_left(c);
		return *this;
	}

	this_type& operator*() {
		return *this;
	}

	this_type& operator++(){
		return *this;
	} //prefix increment
};

}

namespace seedmod {
  	//This is an OutputIterator
    //operator++ <--> ret_m.shift_left(*c);
    //operator=  rembers c
	template <class OutMerIterator,
		      char MATCH_CHR = '#',
		      char DEL_CHR = '^',
			  char INS_CHR = 'v'>
	class SpacedSeedForIndexSquasherIterator
	{
	private:
		char c_;
		const char * seed_;
		OutMerIterator & oit_;

	public:
		typedef SpacedSeedForIndexSquasherIterator<OutMerIterator,
				MATCH_CHR,DEL_CHR,INS_CHR> this_type;

		SpacedSeedForIndexSquasherIterator(const char * seed, OutMerIterator & oit)
			:seed_(seed),oit_(oit){}

		this_type& operator=(const char c) {
			c_=c;
			return *this;
		}

		this_type& operator*() {
			return *this;
		}

		this_type& operator++(){
			do{
				if(DEL_CHR==*seed_){
					++seed_;
					continue;
				}else if(MATCH_CHR==*seed_){
					*oit_ = c_;
				}else if(INS_CHR==*seed_){
				}
				//space
				break;
			}while (true);
			++seed_;
			return *this;
		} //prefix increment


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
	template <class OutMerIterator,
			      char MATCH_CHR = '#',
			      char DEL_CHR = '^',
				  char INS_CHR = 'v'>
	class SpacedSeedForReadSquasherIterator:
			public SpacedSeedForReadSquasherIterator<OutMerIterator,
													MATCH_CHR,INS_CHR,DEL_CHR>
	{
		typedef SpacedSeedForReadSquasherIterator<OutMerIterator,
				MATCH_CHR,DEL_CHR,INS_CHR> this_type;
		typedef SpacedSeedForReadSquasherIterator<OutMerIterator,
				MATCH_CHR,INS_CHR,DEL_CHR> super;

		SpacedSeedForReadSquasherIterator(const char * seed, OutMerIterator & oit)
			:super(seed,oit){}
	};

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

