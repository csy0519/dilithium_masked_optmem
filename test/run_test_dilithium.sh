export PQCLEAN_ONLY_TESTS=common,common_lib,functest,functest_sanitizers,nistkat,testvectors,valgrid
#export PQCLEAN_ONLY_TESTS=functest,functest_sanitizers
export PQCLEAN_ONLY_SCHEMES=dilithium2
export PQCLEAN_ONLY_IMPLEMENTATIONS=mask_optmem

pytest

#python3 test_functest.py


