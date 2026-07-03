function fs_expect(cond, msg)
%FS_EXPECT  Simple assertion: record pass/fail, continue on failure.
%
%   fs_expect(cond, msg)
%       cond — logical; msg — string describing the assertion
%
%   Uses appdata slot 'fs_test_state' to accumulate pass/fail counts.
%   fs_setup initialises it; fs_teardown reports + exits.

    st = getappdata(0, 'fs_test_state');
    if isempty(st)
        st = struct('pass', 0, 'fail', 0, 'first_fail', '');
    end
    if cond
        st.pass = st.pass + 1;
        fprintf('  [ok]   %s\n', msg);
    else
        st.fail = st.fail + 1;
        if isempty(st.first_fail)
            st.first_fail = msg;
        end
        fprintf('  [FAIL] %s\n', msg);
    end
    setappdata(0, 'fs_test_state', st);
end
