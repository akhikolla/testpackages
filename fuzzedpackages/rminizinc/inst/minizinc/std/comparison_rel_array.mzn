%-----------------------------------------------------------------------------%
% Reflect an array of comparison values onto a comparison value variable using
% a lexicographic interpretation of the array.  The comparison values are
% encoded as follows:  > | = | <  as  -1 | 0 | +1.
% Uses of this constraint are generated by Cadmium transformations that
% simplify ordering constraints on expressions of complex types.
%-----------------------------------------------------------------------------%
predicate comparison_rel_array(array[int] of var -1..1: rels, var -1..1: rel) =
          let { var int: firstneg = arg_min(rels);
                var int: firstpos = arg_max(rels); } in
          % if they only consist of one number then return it
          if firstneg = firstpos then
	        rel = rels[min(index_set(rels))]
          else
                %% if there are no -1 then it must be +1
                if rels[firstneg] = 0 then rel = +1 
                %% otherwise if ther are no +1 it must be -1
                elseif rels[firstpos] = 0 then rel = -1 
                %% otherwise if a neg appears earlier than -1
                elseif firstneg < firstpos then rel = -1
                %% otherwise +1
                else rel = +1
          endif;
