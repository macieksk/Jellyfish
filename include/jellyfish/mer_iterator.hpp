/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __MER_ITERATOR_HPP__
#define __MER_ITERATOR_HPP__

#include <iterator>

#include <seedmod/seedmod.hpp>

namespace jellyfish {
template<typename SequencePool, typename MerType, typename FullSeedMerType>
class mer_iterator : public std::iterator<std::input_iterator_tag,MerType> {
  typename SequencePool::job* job_;
  const char*                 cseq_;
  MerType                     ret_m_; // Returned squashed seed mer
  MerType                     ret_rcm_; // Returned squashed seed reverse complement mer
  FullSeedMerType             m_; // mer
  FullSeedMerType             rcm_; // reverse complement mer
  unsigned int                filled_;
  const bool                  canonical_;
  const char *				  seed_;

public:
  typedef MerType      mer_type;
  typedef FullSeedMerType      fullseedmer_type;
  typedef SequencePool sequence_parser_type;
  typedef jellyfish::mer_shift_left_output_iterator<MerType> mer_sleft_oiter;
  typedef seedmod::SpacedSeedForIndexSquasherIterator<char,mer_sleft_oiter> seed_squasher_iter_type;

  mer_iterator(SequencePool& seq, bool canonical = false, const char * seed = NULL) :
    job_(new typename SequencePool::job(seq)), cseq_(0), filled_(0), canonical_(canonical),
	seed_(seed)
  {
    if(job_->is_empty()) {
      delete job_;
      job_ = 0;
    } else {
      cseq_ = (*job_)->start;
      this->operator++();
    }
  }
  mer_iterator() : job_(0), cseq_(0), filled_(0), canonical_(false) { }
  //  mer_iterator(const mer_iterator& rhs) : job_(rhs.job_), cseq_(rhs.cseq_), m_(rhs.m_), filled_(rhs.filled_) { }
  ~mer_iterator() {
    delete job_;
  }

  bool operator==(const mer_iterator& rhs) const { return job_ == rhs.job_; }
  bool operator!=(const mer_iterator& rhs) const { return job_ != rhs.job_; }

  operator void*() const { return (void*)job_; }
  //const mer_type& operator*() const { return !canonical_ || m_ < rcm_ ? m_ : rcm_; }
  const mer_type& operator*() {
	    //Squash full seed mer to squashed mer
	  	auto squash_mer = [&](const fullseedmer_type & fmer,mer_type & ret_m){
	  		mer_sleft_oiter meroiter(ret_m);
	  		seed_squasher_iter_type seed_squash_it(seed_,meroiter);
	  		fmer.to_chars(seed_squash_it);
	  	};
	  	//if (strlen(seed_)>0) {
	  	//if (*seed_!=0) {   //TODO change to template instantiation,
	  						//Now off -- we work only with seeds
	  	squash_mer(m_,ret_m_);
	  	if (canonical_) {
	  		//no reverse complements for asymmetric seeds
	  		//squash_mer(rcm_,ret_rcm_);
	  		std::cerr << std::endl << "Jellyfish-seed does not support canonical k-mers!";
	  		exit(-1);
	  	}

	  	//} else {
	  	//	ret_m_  = m_;
	  	//	ret_rcm_ = rcm_;
	  	//}
	    //std::cerr << std::endl << m_.to_str() << " <--> " << ret_m_.to_str();
	    //return !canonical_ || ret_m_ < ret_rcm_ ? ret_m_ : ret_rcm_;
	  	return ret_m_;
  }
  //const mer_type& fullmer() const { return !canonical_ || ret_m_ < ret_rcm_ ? m_ : rcm_; }
  //bool canonical_is_reversed() const { return !canonical_ || ret_m_ < ret_rcm_ ? false : true; }
  const mer_type* operator->() { return &this->operator*(); }
  mer_iterator& operator++() {
    while(true) {
      while(cseq_ == (*job_)->end) {
        job_->next();
        if(job_->is_empty()) {
          delete job_;
          job_  = 0;
          cseq_ = 0;
          return *this;
        }
        cseq_   = (*job_)->start;
        filled_ = 0;
      }

      do {
        int code = m_.code(*cseq_++);
        if(code >= 0) {
          m_.shift_left(code);
          if(canonical_)
            rcm_.shift_right(rcm_.complement(code));
          filled_ = std::min(filled_ + 1, fullseedmer_type::k());
        } else
          filled_ = 0;
      } while(filled_ < m_.k() && cseq_ < (*job_)->end);
      if(filled_ >= m_.k())
        break;
    }


    return *this;
  }

  mer_iterator operator++(int) {
    mer_iterator res(*this);
    ++*this;
    return res;
  }
};

} // namespace jellyfish {

#endif /* __MER_ITERATOR_HPP__ */
