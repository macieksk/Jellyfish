/*
 * seedmod.hpp
 *
 *  Created on: Dec 16, 2014
 *      Author: Maciek Sykulski <macieksk@gmail.com>
 */

#ifndef SEEDMOD_HPP_
#define SEEDMOD_HPP_

#include <string.h>

namespace jellyfish {
  namespace seedmod {


  	//This is an OutputIterator
    //operator++ <--> ret_m.shift_left(*c);
    //operator=  rembers c
	template <class MerTypeOut,
			 char MATCH_CHR = '#'>
	class SpacedSeedSquasherIterator : public std::vector<char>
	{
		//static const char match_chr = MATCH_CHR;
		typedef SpacedSeedSquasherIterator<MerTypeOut,MATCH_CHR> this_type;

	public:
		SpacedSeedSquasherIterator(const char * seed, MerTypeOut & mer)
		:seed_(seed),mer_(mer){}

		this_type& operator=(const char c) {
			c_=c;
			return *this;
		}

		this_type& operator*() {
			return *this;
		}

		this_type& operator++(){
			if (*seed_==MATCH_CHR)
				mer_.shift_left(c_);
			++seed_;
			return *this;
		} //prefix increment


		bool at_end() const{
			return seed_==0;
		}


	private:
		char c_;
		const char * seed_;
		MerTypeOut & mer_;

	public:

		static size_t weight(const char * s){
			size_t w=0;
			for (;*s;++s)
				if (*s==MATCH_CHR) ++w;
			return w;
		}
		static size_t span(const char * s){
			return strlen(s);
		}

	};


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

