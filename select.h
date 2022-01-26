
#ifndef _IBPEP_SELECT_H_
#define _IBPEP_SELECT_H_

#include "char_io.h"
#include "pattern_match.h"

PatternMatch *pick_pattern_index(CharacterProducer * const & ff,
				 int pmselect,
				 int nmismatch,
				 std::vector<std::pair<int,int> > *exact_const,
				 std::vector<int> *patlen,
				 int seedlen,
				 bool wildcard,
				 bool textn,
				 bool indels,
				 bool dna_mut,
				 char eos,
				 bool verbose);
#endif
