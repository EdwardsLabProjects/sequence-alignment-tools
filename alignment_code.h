
#ifndef _IBPEP_alignment_code_h
#define _IBPEP_alignment_code_h

# define alignment_codes 13
enum alignment_code {
  alignment_none=0,
  alignment_equal=1,
  alignment_wildcard_equal=2,
  alignment_substitution=3,
  alignment_insertion=4,
  alignment_deletion=5,
  alignment_constraint_violation=6,
  alignment_end = 7,
  alignment_substitution_1=8,
  alignment_substitution_2=9,
  alignment_substitution_3=10,
  alignment_insertion_3=11,
  alignment_deletion_3=12
};

enum alignment_masks {
  alignment_mask_none = 1,
  alignment_mask_equal = 2,
  alignment_mask_wildcard_equal = 4,
  alignment_mask_substitution = 8,
  alignment_mask_insertion = 16,
  alignment_mask_deletion = 32,
  alignment_mask_constraint_violation = 64,
  alignment_mask_alignment_end = 128,
  alignment_mask_substitution_1 = 256,
  alignment_mask_substitution_2 = 512,
  alignment_mask_substitution_3 = 1024,
  alignment_mask_insertion_3 = 2048,
  alignment_mask_deletion_3 = 4096
};

#endif
