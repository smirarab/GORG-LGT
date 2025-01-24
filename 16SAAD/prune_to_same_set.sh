python ~/my_gits/myTools/prune_tree.py -i astral.tre -l "`nw_labels -I aad_tree.nwk`" -v -o astral_pruned.tre
sed -i -e "s/'//g" astral_pruned.tre
python ~/my_gits/myTools/prune_tree.py -i aad_tree.nwk -l "`nw_labels -I astral_pruned.tre`" -v -o aad_tree_pruned.nwk
sed -i -e "s/'//g" aad_tree_pruned.nwk
python ~/my_gits/myTools/prune_tree.py -i RAxML.T16S-SINA-1000-bootcutoff_MVrooted -l "`nw_labels -I astral_pruned.tre`" -v -o RAxML.T16S-SINA-1000-bootcutoff_MVrooted_pruned
sed -i -e "s/'//g" RAxML.T16S-SINA-1000-bootcutoff_MVrooted_pruned
python ~/my_gits/myTools/prune_tree.py -i RAxML.T23S-SINA-1000-bootcutoff_MVrooted -l "`nw_labels -I astral_pruned.tre`" -v -o RAxML.T23S-SINA-1000-bootcutoff_MVrooted_pruned
sed -i -e "s/'//g" RAxML.T23S-SINA-1000-bootcutoff_MVrooted_pruned
