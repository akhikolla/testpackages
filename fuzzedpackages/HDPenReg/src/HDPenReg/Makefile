all:
	($(MAKE) -f lassoModels.mk)
	($(MAKE) -f lars.mk)

clean:
	($(RM) lassoModels/*.o)
	($(RM) ./*.o)
	($(RM) lars/*.o)

# use to debug
echo:
	@echo $(ALL_CXXFLAGS)
	@echo $(ALL_CPPFLAGS)
	@echo $(HD_CXXFLAGS)
	@echo $(HD_CPPFLAGS)
